from Bio.SeqUtils.MeltingTemp import Tm_staluc
from primer3.bindings import calcHairpin, calcHomodimer, calcHeterodimer
import re
from collections import defaultdict
from util import grouper

dG_SCALE = 1e3
complementary_bases = {"a": "t", "t": "a", "c": "g", "g": "c"}

# FROM: http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333&aa=1&style=N
# fields: [triplet] [amino acid] [fraction] [frequency: per thousand] ([number])
codon_usage = """
UUU F 0.57 19.7 (   101)  UCU S 0.11  5.7 (    29)  UAU Y 0.53 16.8 (    86)  UGU C 0.42  5.9 (    30)
UUC F 0.43 15.0 (    77)  UCC S 0.11  5.5 (    28)  UAC Y 0.47 14.6 (    75)  UGC C 0.58  8.0 (    41)
UUA L 0.15 15.2 (    78)  UCA S 0.15  7.8 (    40)  UAA * 0.64  1.8 (     9)  UGA * 0.36  1.0 (     5)
UUG L 0.12 11.9 (    61)  UCG S 0.16  8.0 (    41)  UAG * 0.00  0.0 (     0)  UGG W 1.00 10.7 (    55)

CUU L 0.12 11.9 (    61)  CCU P 0.17  8.4 (    43)  CAU H 0.55 15.8 (    81)  CGU R 0.36 21.1 (   108)
CUC L 0.10 10.5 (    54)  CCC P 0.13  6.4 (    33)  CAC H 0.45 13.1 (    67)  CGC R 0.44 26.0 (   133)
CUA L 0.05  5.3 (    27)  CCA P 0.14  6.6 (    34)  CAA Q 0.30 12.1 (    62)  CGA R 0.07  4.3 (    22)
CUG L 0.46 46.9 (   240)  CCG P 0.55 26.7 (   137)  CAG Q 0.70 27.7 (   142)  CGG R 0.07  4.1 (    21)

AUU I 0.58 30.5 (   156)  ACU T 0.16  8.0 (    41)  AAU N 0.47 21.9 (   112)  AGU S 0.14  7.2 (    37)
AUC I 0.35 18.2 (    93)  ACC T 0.47 22.8 (   117)  AAC N 0.53 24.4 (   125)  AGC S 0.33 16.6 (    85)
AUA I 0.07  3.7 (    19)  ACA T 0.13  6.4 (    33)  AAA K 0.73 33.2 (   170)  AGA R 0.02  1.4 (     7)
AUG M 1.00 24.8 (   127)  ACG T 0.24 11.5 (    59)  AAG K 0.27 12.1 (    62)  AGG R 0.03  1.6 (     8)

GUU V 0.25 16.8 (    86)  GCU A 0.11 10.7 (    55)  GAU D 0.65 37.9 (   194)  GGU G 0.29 21.3 (   109)
GUC V 0.18 11.7 (    60)  GCC A 0.31 31.6 (   162)  GAC D 0.35 20.5 (   105)  GGC G 0.46 33.4 (   171)
GUA V 0.17 11.5 (    59)  GCA A 0.21 21.1 (   108)  GAA E 0.70 43.7 (   224)  GGA G 0.13  9.2 (    47)
GUG V 0.40 26.4 (   135)  GCG A 0.38 38.5 (   197)  GAG E 0.30 18.4 (    94)  GGG G 0.12  8.6 (    44)
"""

codon_to_aa = {}
aa_to_codon = defaultdict(list)
codon_count = {}

for codon, aa, _, _, count in grouper(re.sub("\(|\)", "", codon_usage).split(), 5):
    codon = codon.lower().replace("u", "t")
    codon_to_aa[codon] = aa
    aa_to_codon[aa].append(codon)
    codon_count[codon] = int(count)
aa_to_codon = dict(aa_to_codon)

codons_by_usage = sorted(codon_count.keys(), key=codon_count.__getitem__, reverse=True)

for aa, codons in aa_to_codon.items():
    aa_to_codon[aa] = sorted(codons, key=codon_count.__getitem__, reverse=True)


def site_diff(a, b):
    diff = []
    for i, (x_a, x_b) in enumerate(zip(a, b)):
        if x_a != x_b:
            diff.append((i, x_a, x_b))
    return diff


def design_tm(seq, start, tm_goal, reverse=False, gc_clamp=True):
    stop = start + 1
    while 0 <= stop <= len(seq) and 0 <= start < len(seq):
        subseq = seq[start:stop]
        tm = Tm_staluc(subseq)
        tm_goal_met = False
        gc_clamp_goal_met = False
        if not gc_clamp:
            gc_clamp_goal_met = True
        else:
            if reverse and subseq[:1] in ("g", "c"):
                gc_clamp_goal_met = True
            elif not reverse and subseq[-1] in ("g", "c"):
                gc_clamp_goal_met = True
        if isinstance(tm_goal, tuple):
            if tm_goal[0] <= tm <= tm_goal[1]:
                tm_goal_met = True
        else:
            if tm_goal <= tm:
                tm_goal_met = True
        if tm_goal_met and gc_clamp_goal_met:
            return subseq, tm
        if reverse:
            start -= 1
        else:
            stop += 1
    raise ValueError("ran into end of sequence before achieved goal")


def reverse_complement(seq):
    return "".join(complementary_bases.get(base, base) for base in reversed(seq))


def iva_substitution_primers(
    seq,
    mutation,
    start,
    stop=None,
    binding_tm=(59, 62),
    homology_tm=(48, 52),
    thermodynamics=True,
):
    if stop is None:
        stop = start + len(mutation)
    reversion_mutation = seq[start:stop]
    try:
        forward_binding, forward_binding_tm = design_tm(
            seq, stop, binding_tm, gc_clamp=True
        )
        forward_homology, forward_homology_tm = design_tm(
            seq, start - 1, homology_tm, reverse=True
        )
        reverse_primer, reverse_tm = design_tm(
            seq, start - 1, binding_tm, reverse=True, gc_clamp=True
        )
    except:
        return {}
    forward_primer = forward_homology + mutation + forward_binding
    reversion_forward_primer = forward_homology + reversion_mutation + forward_binding
    reverse_primer = reverse_complement(reverse_primer)
    res = {
        "forward_primer": forward_primer,
        "reverse_primer": reverse_primer,
        "reversion_forward_primer": reversion_forward_primer,
        "forward_homology_tm": forward_homology_tm,
        "forward_binding_tm": forward_binding_tm,
        "reverse_tm": reverse_tm,
        "reverse_overhang": len(reverse_primer) - len(forward_homology),
        "forward_homology_len": len(forward_homology),
        "forward_binding_len": len(forward_binding),
    }
    if thermodynamics:
        for seq, prefix in ((forward_primer, "forward_"), (reverse_primer, "reverse_")):
            if len(seq) > 60:
                continue
            homodimer_res = calcHomodimer(seq)
            res[prefix + "homodimer_tm"] = homodimer_res.tm
            res[prefix + "homodimer_dG"] = homodimer_res.dg / dG_SCALE
            hairpin_res = calcHairpin(seq)
            res[prefix + "hairpin_tm"] = hairpin_res.tm
            res[prefix + "hairpin_dG"] = hairpin_res.dg / dG_SCALE
        heterodimer_res = calcHeterodimer(forward_primer, reverse_primer)
        res["heterodimer_tm"] = heterodimer_res.tm
        res["heterodimer_dG"] = heterodimer_res.dg / dG_SCALE
    return res


def mutation_name(old, new, position):
    if codon_to_aa[old] != "*":
        old = codon_to_aa[old].upper()
    if codon_to_aa[new] != "*":
        new = codon_to_aa[new].upper()
    return "{}{:d}{}".format(old, position, new)


def parse_mutation_name(s):
    from_aa, res, to_aa = re.match("([A-Zatcg]+)(\d+)([A-Zatcg]+)", s).groups()
    res = int(res)
    return from_aa, to_aa, res


def find_snps(seq, mutations, synonymous=False):
    snps = []
    for (res, from_codon, to_codon) in mutations:
        res0 = res - 1  # 0-indexed residue number
        if synonymous:
            from_aa = codon_to_aa[seq[res0 : res0 + 3].lower()]
            alt_from_codons = aa_to_codon[from_aa]
        else:
            from_codons = [from_codon]
        for alt_from_codon in from_codons:
            diff = site_diff(alt_from_codon, to_codon)
            if len(diff) == 1:
                snp = {
                    "residue": res,
                    "from_codon": from_codon,
                    "to_codon": to_codon,
                    "aa_mutation_name": mutation_name(alt_from_codon, to_codon, res),
                    "aa_mutation": (res, alt_from_codon, to_codon),
                    "nt_mutation": (3 * res0 + diff[0][0], diff[0][1], diff[0][2]),
                    "synonymous": alt_from_codon != from_codon,
                    "transition": "{}->{}".format(diff[0][1], diff[0][2]),
                }
                snps.append(snp)
    return snps


def stop_codon_snps(seq, synonymous=False):
    stop_codons = aa_to_codon["*"]
    mutations = []
    for res0, from_codon in enumerate(grouper(seq, 3)):
        from_codon = "".join(from_codon)
        res = res0 + 1
        for stop_codon in stop_codons:
            mutations.append((res, from_codon, stop_codon))
    return find_snps(seq, mutations, synonymous=synonymous)


def stop_codon_primers(seq, synonymous=False):
    stop_codons = aa_to_codon["*"]
    primers = []
    for snp in stop_codon_snps(seq, synonymous=synonymous):
        primer = iva_substitution_primers(
            seq, snp["nt_mutation"][2], seq["nt_mutation"][0]
        )
        primer.update(snp)
        primers.append(primer)
    return primers
