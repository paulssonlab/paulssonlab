from cytoolz import partial
from lxml import etree
from requests_ratelimiter import LimiterSession
from paulssonlab.cloning.io import bytes_to_value
from paulssonlab.util import first, only

ENTREZ_BASE_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# aliases for commonly-used organisms
ORGANISMS = {
    "ecoli": "Escherichia coli str. K-12 substr. MG1655",
    "bsubtilis": "Bacillus subtilis subsp. subtilis str. 168",
}

session = LimiterSession(per_second=3)


def escape(s):
    return s.replace('"', '\\"')


def get(util, method="get", retmode="xml", **params):
    url = f"{ENTREZ_BASE_URL}/{util}.fcgi"
    params = {"retmode": retmode, **params}
    if method == "get":
        res = session.get(url, params=params)
    elif method == "post":
        res = session.post(url, params=params)
    else:
        raise ValueError(f"unknown method '{method}'")
    res.raise_for_status()
    if retmode == "xml":
        # print("***")
        # print(res.content.decode())
        # print("***")
        return etree.fromstring(res.content)
    else:
        return res.content


def get_gene(name, organism, one=True):
    organism = ORGANISMS.get(organism, organism)
    # we restrict gene results to those with RefSeq entries
    term = f'{escape(name)}[Gene/Protein Name] AND "{escape(organism)}"[Organism] AND "srcdb refseq"[Properties] AND alive[prop]'
    search_results = get(
        "esearch",
        db="gene",
        term=term,
    )
    ids = search_results.xpath("./IdList/Id/text()")
    if len(ids) == 0:
        raise ValueError(f"could not find '{name}' (organism: '{organism}')")
    if one:
        if len(ids) != 1:
            raise ValueError(f"found {len(ids)} ids, expecting 1")
    entries = get("efetch", db="gene", id=",".join(ids))
    if one:
        return entries[0]
    else:
        return entries


def get_gene_seq_loc(name, organism, rettype="fasta", one=True):
    genes = get_gene(name, organism, one=one)
    if one:
        genes = [genes]
    locs = []
    for gene in genes:
        commentary = only(
            (
                l
                for l in gene.findall("Entrezgene_locus/Gene-commentary")
                if l.findtext("Gene-commentary_label") != "old_locus_tag"
            ),
            msg="expecting a single 'Gene-commentary' not marked with old_locus_tag",
        )
        # only handle case with a single genome interval
        interval = only(
            commentary.findall(
                "Gene-commentary_products/Gene-commentary/Gene-commentary_genomic-coords/Seq-loc/Seq-loc_int/Seq-interval"
            ),
            msg="expecting a single 'Seq-interval'",
        )
        gi = int(interval.findtext("Seq-interval_id/Seq-id/Seq-id_gi"))
        interval_from = int(interval.findtext("Seq-interval_from"))
        interval_to = int(interval.findtext("Seq-interval_to"))
        locs.append((gi, interval_from, interval_to))
    return locs


def get_gene_seq(name, organism, rettype="fasta", one=True):
    locs = get_gene_seq_loc(name, organism, rettype=rettype, one=one)
    seqs = []
    for gi, interval_from, interval_to in locs:
        # interval_from/to here seem to be 0-indexed,
        # need to increment by 1 before passing to efetch
        seq_raw = get(
            "efetch",
            db="nuccore",
            id=gi,
            rettype=rettype,
            retmode="text",
            **{"from": interval_from + 1, "to": interval_to + 1},
        )
        # SEE: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
        if rettype in ("fasta", "fasta_cds_na", "fasta_cds_aa"):
            seq = bytes_to_value(seq_raw, "chemical/seq-na-fasta")
        elif rettype in ("gb", "gbwithparts"):
            seq = bytes_to_value(seq_raw, "chemical/seq-na-genbank")
        else:
            raise ValueError(f"unknown rettype '{rettype}'")
        seqs.append(seq)
    if one:
        return seqs[0]
    else:
        return seqs
