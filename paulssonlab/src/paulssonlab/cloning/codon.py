import pandas as pd
import requests
import re
from operator import itemgetter
from paulssonlab.util import grouper, first

KAZUSA_URI = "http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?"


def get_codon_usage_table(species_id=37762, replace_uracil=True):
    """Get a codon usage table from the Kazusa web service.

    Defaults to E. coli.
    """
    table = _get_codon_usage_table(species_id)
    return parse_codon_usage_table(table)


def _get_codon_usage_table(species_id):
    r = requests.get(KAZUSA_URI, params={"species": species_id, "aa": 1, "style": "N"})
    r.raise_for_status()
    match = re.search(r"<PRE>(.*)</PRE>", r.content.decode(), re.DOTALL)
    if not match:
        raise Exception("could not find codon table in webpage")
    return match.group(1)


def parse_codon_usage_table(table, replace_uracil=True):
    rows = []
    for codon, aa, _, _, count in grouper(re.sub("\(|\)", "", table).split(), 5):
        if replace_uracil:
            codon = codon.replace("U", "T")
        count = int(count)
        rows.append((codon, aa, count))
    # return pd.DataFrame(rows, columns=["codon", "aa", "count"])
    return rows


def aa_frequency(codon_usage_table=None):
    if codon_usage_table is None:
        codon_usage_table = get_codon_usage_table()
    aa_freq = {}
    for codon, aa, freq in codon_usage_table:
        if aa not in aa_freq:
            aa_freq[aa] = 0
        aa_freq[aa] += freq
    return aa_freq


def codons_by_relative_frequency(codon_usage_table=None):
    if codon_usage_table is None:
        codon_usage_table = get_codon_usage_table()
    total_counts = sum(t[2] for t in codon_usage_table)
    aa_freq = aa_frequency(codon_usage_table)
    codon_usage_table = sorted(codon_usage_table, key=itemgetter(2), reverse=True)
    aa_to_codons = {}
    for codon, aa, freq in codon_usage_table:
        if aa not in aa_to_codons:
            aa_to_codons[aa] = {}
        aa_to_codons[aa][codon] = freq / aa_freq[aa]
    return aa_to_codons


def codons_by_absolute_frequency(codon_usage_table=None):
    if codon_usage_table is None:
        codon_usage_table = get_codon_usage_table()
    total_counts = sum(t[2] for t in codon_usage_table)
    codon_usage_table = sorted(codon_usage_table, key=itemgetter(2), reverse=True)
    aa_to_codons = {}
    for codon, aa, freq in codon_usage_table:
        if aa not in aa_to_codons:
            aa_to_codons[aa] = {}
        aa_to_codons[aa][codon] = freq / total_counts
    return aa_to_codons


def back_translate(aa_seq, aa_to_codons=None):
    if aa_to_codons is None:
        aa_to_codons = most_common_codons(codon_usage_table)
    seq = "".join([first(aa_to_codons[aa].keys()) for aa in aa_seq])
    return seq
