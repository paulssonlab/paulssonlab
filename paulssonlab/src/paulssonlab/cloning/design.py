import random
from paulssonlab.cloning.sequence import reverse_complement
from paulssonlab.cloning.workflow import normalize_seq


def random_bases(n, letters="atcg", seed=""):
    seed = normalize_seq(seed)
    myrandom = random.Random(seed)
    return "".join(myrandom.choice(letters) for _ in range(max(n, 0)))


def type2s_with_spacer(enzyme, overhang_length):
    return _type2s_with_spacer(enzyme.site, enzyme.fst5, enzyme.fst3, overhang_length)


def _type2s_with_spacer(site, cut5, cut3, overhang_length):
    cut5 = cut5 - len(site)
    if cut5 < 0 or cut3 < 0:
        raise NotImplementedError(
            "cannot handle cuts that intersect with enzyme recognition site"
        )
    if cut5 > cut3:
        # overhang = sequence.reverse_complement(overhang)
        raise NotImplementedError(
            "cannot handle enzymes where 5' cut is downstream of 3' cut"
        )
    if cut3 - cut5 != overhang_length:
        raise ValueError(
            f"expecting an overhang of length {cut3 - cut5}, instead got {overhang_length}"
        )
    seq = site + random_bases(cut5, seed=site)
    return seq


def golden_gate_placeholder(
    inner_enzyme,
    outer_enzyme,
    overhang1,
    overhang2,
    num_random_bases=6,
    random_flanks=True,
):
    seq = overhang1 + reverse_complement(
        type2s_with_spacer(inner_enzyme, len(overhang1))
    )
    seq += (
        random_bases(2 * num_random_bases, seed=seq)
        + type2s_with_spacer(inner_enzyme, len(overhang2))
        + overhang2
    )
    if outer_enzyme is not None:
        seq = (
            type2s_with_spacer(outer_enzyme, len(overhang1))
            + seq
            + reverse_complement(type2s_with_spacer(outer_enzyme, len(overhang2)))
        )
        if random_flanks:
            seq = random_bases(num_random_bases, seed=seq) + seq
            seq = seq + random_bases(num_random_bases, seed=seq)
    return seq
