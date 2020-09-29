import re
from datetime import datetime, date, timezone
from copy import deepcopy
from Bio.SeqFeature import SeqFeature, FeatureLocation
from paulssonlab.cloning.sequence import (
    slice_seq,
    reverse_complement,
    join_seqs,
    get_seq,
)
from paulssonlab.util import sign


def _re_search(enzyme, seq, linear=True):
    compsite = re.compile(
        enzyme.compsite.pattern, enzyme.compsite.flags | re.IGNORECASE
    )
    if not linear:
        seq = seq + seq[1 : enzyme.size]
    re_sites = [
        (i.start(), i.group(1) is not None)
        for i in re.finditer(compsite, str(get_seq(seq)))
    ]
    return re_sites


def _re_search_cuts(binding_locs, enzyme):
    cuts = []
    for loc, cut_upstream in binding_locs:
        for cut5, cut3 in ((enzyme.fst5, enzyme.fst3), (enzyme.scd5, enzyme.scd3)):
            if cut5 is None and cut3 is None:
                continue
            if cut_upstream:
                if cut5 is not None:
                    cut5_loc = loc + cut5
                else:
                    cut5_loc = None
                if cut3 is not None:
                    cut3_loc = loc + enzyme.size + cut3
                else:
                    cut3_loc = None
            else:
                if cut3 is not None:
                    cut5_loc = loc - cut3
                else:
                    cut5_loc = None
                if cut5 is not None:
                    cut3_loc = loc + enzyme.size - cut5
                else:
                    cut3_loc = None
            # is_5prime_overhang is true if cut5 is upstream of cut3
            if cut5 is not None and cut3 is not None:
                # is_5prime_overhang = cut5_loc > cut3_loc
                sense = sign(cut3_loc - cut5_loc)
            else:
                # is_5prime_overhang = None
                sense = 0
            cuts.append((cut5_loc, cut3_loc, sense, cut_upstream))
    return cuts


def re_search(seq, enzyme, linear=True):
    binding_locs = _re_search(enzyme, seq, linear=linear)
    cuts = _re_search_cuts(binding_locs, enzyme)
    length = len(seq)
    return sorted([(c[0] % length, c[1] % length, *c[2:]) for c in cuts])


def _get_overhang(seq, cut5, cut3, sense, cut_upstream):
    if (sense == -1 and not cut_upstream) or (sense == 1 and cut_upstream):
        loc = cut3
    else:
        loc = cut5
    seq = get_seq(seq)  # get Seq from SeqRecord, if necessary
    return ((slice_seq(seq, cut5, cut3), sense), loc)


def _re_digest(seq, cuts):
    cuts.append(cuts[0])
    seqs = []
    for cut1, cut2 in zip(cuts[:-1], cuts[1:]):
        # check for sequences with inward-facing RE binding sites
        if cut1[3] == True and cut2[3] == False:
            overhang1, loc1 = _get_overhang(seq, *cut1)
            overhang2, loc2 = _get_overhang(seq, *cut2)
            seq = slice_seq(seq, loc1, loc2)
            seqs.append((seq, overhang1, overhang2))
    seqs = sorted(seqs, key=lambda x: len(x[0]))
    return seqs


def re_digest(seq, enzyme, linear=False):
    cuts = re_search(seq, enzyme, linear=linear)
    return _re_digest(seq, cuts)


def _check_seq_compatibility(seq1, seq2):
    _, overhang1_1, overhang1_2 = seq1
    _, overhang2_1, overhang2_2 = seq2
    return (overhang1_2[0].upper(), overhang1_2[1]) == (
        overhang2_1[0].upper(),
        overhang2_1[1],
    )


def _reverse_complement_overhangs(seq_with_overhangs):
    seq, overhang1, overhang2 = seq_with_overhangs
    overhang1_rc = (reverse_complement(overhang1[0]), overhang1[1])
    overhang2_rc = (reverse_complement(overhang2[0]), overhang2[1])
    return (reverse_complement(seq), overhang2_rc, overhang1_rc)


def _5prime_overhang(overhang):
    if not overhang[1]:
        return reverse_complement(overhang[0])
    else:
        return overhang[0]


def ligate(seqs, linear=True):
    if len(seqs) < 2:
        raise ValueError("need at least two sequences to assemble")
    seq1 = seqs[0]
    seq2 = seqs[1]
    seq1_rc = _reverse_complement_overhangs(seq1)
    seq2_rc = _reverse_complement_overhangs(seq2)
    if _check_seq_compatibility(seq1, seq2):
        pass
    elif _check_seq_compatibility(seq1, seq2_rc):
        seqs[1] = seq2_rc
    elif _check_seq_compatibility(seq1_rc, seq2):
        seqs[0] = seq1_rc
    elif _check_seq_compatibility(seq1_rc, seq2_rc):
        seqs[0] = seq1_rc
        seqs[1] = seq2_rc
    else:
        raise ValueError(f"overhang mismatch when assembling sequences 0 and 1")
    seqs = [*seqs, None]
    seqs_to_join = []
    for idx in range(len(seqs) - 1):
        seq1 = seqs[idx]
        seq2 = seqs[idx + 1]
        if seq2 is not None:
            seq2_rc = _reverse_complement_overhangs(seq2)
            if _check_seq_compatibility(seq1, seq2):
                pass
            elif _check_seq_compatibility(seq1, seq2_rc):
                # raise NotImplementedError
                seq2 = seqs[idx + 1] = seq2_rc  # TODO: does this change zip?
            else:
                raise ValueError(
                    f"overhang mismatch when assembling sequences {idx} and {idx + 1}: {seq1[2]} does not match {seq2[1]} or {seq2_rc[1]}"
                )
        seqs_to_join.append(_5prime_overhang(seq1[1]))
        seqs_to_join.append(seq1[0])
        if seq2 is None:
            if linear:
                seqs_to_join.append(_5prime_overhang(seq1[2]))
            else:
                if not _check_seq_compatibility(seq1, seqs[0]):
                    raise ValueError(
                        f"overhang mismatch when assembling sequences {idx} and 0 (circularizing): {seq1[2]} does not match {seqs[0][1]}"
                    )
    # copy SeqRecords, add annotations for each part?? (including overhangs)
    joined_seq = join_seqs(seqs_to_join)
    joined_seq.annotations = {
        "molecule_type": "ds-DNA",
        "data_file_division": "SYN",
        "date": datetime.now(timezone.utc).isoformat(),
        "source": "synthetic DNA construct",
        "organism": "synthetic DNA construct",
    }
    if linear is False:
        joined_seq.annotations["topology"] = "circular"
    elif linear:
        joined_seq.annotations["topology"] = "linear"
    # add circular annotation?
    return joined_seq


def assemble(seqs, linear=True):
    seqs_to_assemble = []
    for seq, enzyme, part_name, part_type in seqs:
        subseqs = re_digest(seq, enzyme, linear=linear)
        part, overhang1, overhang2 = subseqs[0]
        if hasattr(part, "features"):
            part = deepcopy(part)
            label = SeqFeature(
                FeatureLocation(-len(overhang1[0]), len(part) + len(overhang2[0])),
                type=part_type,
            )
            label.qualifiers["label"] = [part_name]
            features = [
                feature for feature in part.features if feature.type != "source"
            ]
            part.features = [label, *features]
        seqs_to_assemble.append((part, overhang1, overhang2))
    assembly = ligate(seqs_to_assemble, linear=linear)
    return assembly
