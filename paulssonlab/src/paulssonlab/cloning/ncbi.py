import Bio.Entrez as Entrez
from cytoolz import partial
from paulssonlab.cloning.io import bytes_to_value
from paulssonlab.util import first, only

Entrez.email = "paulssonlab@gmail.com"


def _escape(s):
    return s.replace('"', '\\"')


def get_gene(name, organism, one=True):
    # we restrict gene results to those with RefSeq entries
    search_results = Entrez.read(
        Entrez.esearch(
            "gene",
            f'{_escape(name)}[Gene/Protein Name] AND "{_escape(organism)}"[Organism] AND "srcdb refseq"[Properties] AND alive[prop]',
        )
    )
    ids = search_results.get("IdList", [])
    if len(ids) == 0:
        raise ValueError(f"could not find '{name}' (organism: '{organism}')")
    if one:
        if len(ids) != 1:
            raise ValueError(f"found {len(ids)} ids, expecting 1")
    entries = list(Entrez.parse(Entrez.efetch("gene", id=",".join(ids), retmode="xml")))
    if one:
        return entries[0]
    else:
        return entries


def get_gene_seq(name, organism, rettype="fasta", one=True):
    genes = get_gene(name, organism, one=one)
    if one:
        genes = [genes]
    seqs = []
    for gene in genes:
        locus = only(
            l
            for l in gene["Entrezgene_locus"]
            if l.get("Gene-commentary_label") != "old_locus_tag"
        )
        # only handle case with a single genome interval
        interval = only(
            only(
                locus["Gene-commentary_products"],
                msg="expecting a single 'Gene-commentary_products'",
            )["Gene-commentary_genomic-coords"],
            msg="expecting a single 'Seq-loc_int'",
        )["Seq-loc_int"]["Seq-interval"]
        interval_from = int(interval["Seq-interval_from"])
        interval_to = int(interval["Seq-interval_to"])
        gi = interval["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
        return interval_from, interval_to, gi
        # interval_from/to here seem to be 0-indexed,
        # need to increment by 1 before passing to efetch
        seq_raw = Entrez.efetch(
            db="nuccore",
            id=gi,
            rettype=rettype,
            **{"from": interval_from + 1, "to": interval_to + 1},
        ).read()
        if rettype == "fasta":
            seq = bytes_to_value(seq_raw, "chemical/seq-na-fasta")
        elif rettype == "gb":
            seq = bytes_to_value(seq_raw, "chemical/seq-na-genbank")
        else:
            raise ValueError(f"unknown rettype '{rettype}'")
        seqs.append(seq)
    if one:
        return seqs[0]
    else:
        return seqs
