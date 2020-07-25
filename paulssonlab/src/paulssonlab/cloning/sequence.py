import re
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition
from paulssonlab.util import sign


def reverse_complement(seq):
    if hasattr(seq, "features"):
        return seq.reverse_complement(
            id=True, name=True, description=True, annotations=True, dbxrefs=True
        )
    else:
        return seq.reverse_complement()


def slice_seq(seq, start, end, annotation_start=None, annotation_end=None):
    if start is None:
        start = 0
    if end is None:
        end = len(seq)
    if start < 0:
        start = start % len(seq)
    if end < 0:
        end = end % len(seq)
    if annotation_start is None:
        annotation_start = start
    if annotation_end is None:
        annotation_end = end
    if end < start:
        slice1 = slice_seq(seq, start, None)
        slice2 = slice_seq(seq, None, end)
        new_seq = slice1 + slice2
        # TODO: join features that wrap around? using a join_features func
    else:
        new_seq = seq[start:end]
        if hasattr(seq, "features"):
            features = []
            for feature in seq.features:
                if (
                    annotation_end <= feature.location.nofuzzy_start
                    or annotation_start >= feature.location.nofuzzy_end
                ):
                    continue
                new_feature = feature._shift(-start)
                start_loc = new_feature.location.start
                end_loc = new_feature.location.end
                if annotation_start > feature.location.nofuzzy_start:
                    start_loc = ExactPosition(start - annotation_start)
                if annotation_end < feature.location.nofuzzy_end:
                    end_loc = ExactPosition(annotation_end - start)
                new_feature.location = FeatureLocation(
                    start_loc, end_loc, strand=new_feature.location.strand
                )
                features.append(new_feature)
            new_seq.features = features
        if hasattr(seq, "letter_annotations"):
            # annotation_start and annotation_end only affect features
            # letter_annotations are dropped outside of start:end
            for key, value in seq.letter_annotations.items():
                new_seq._per_letter_annotations[key] = value[start:end]
    return new_seq


def _search_re(enzyme, seq, linear=True):
    compsite = re.compile(
        enzyme.compsite.pattern, enzyme.compsite.flags | re.IGNORECASE
    )
    if not linear:
        seq = seq + seq[1 : enzyme.size]
    re_sites = [
        (i.start(), i.group(1) is not None) for i in re.finditer(compsite, str(seq.seq))
    ]
    return re_sites


def _re_digest_cuts(binding_locs, enzyme):
    cuts = []
    for loc, is_upstream in binding_locs:
        for cut5, cut3 in ((enzyme.fst5, enzyme.fst3), (enzyme.scd5, enzyme.scd3)):
            if cut5 is None and cut3 is None:
                continue
            if is_upstream:
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
            cuts.append((cut5_loc, cut3_loc, sense, is_upstream))
    return cuts


def re_digest(seq, enzyme, linear=True):
    binding_locs = _search_re(enzyme, seq, linear=linear)
    cuts = _re_digest_cuts(binding_locs, enzyme)
    length = len(seq)
    return sorted([(c[0] % length, c[1] % length, *c[2:]) for c in cuts])


def _get_overhang(seq, cut5, cut3, sense):
    if not sense:
        loc = cut5
    else:
        loc = cut3
    if hasattr(seq, "seq"):
        seq = seq.seq  # get Seq if we have a SeqRecord
    return ((slice_seq(seq, cut5, cut3), sense), loc)


def _digest_for_assembly(seq, cuts):
    cuts.append(cuts[0])
    seqs = []
    for cut1, cut2 in zip(cuts[:-1], cuts[1:]):
        # check for sequences with inward-facing RE binding sites
        if cut1[3] == True and cut2[3] == False:
            overhang1, loc1 = _get_overhang(seq, *cut1[:3])
            overhang2, loc2 = _get_overhang(seq, *cut2[:3])
            seq = slice_seq(seq, loc1, loc2)
            seqs.append((seq, overhang1, overhang2))
    seqs = sorted(seqs, key=lambda x: len(x[0]))
    return seqs


def digest_for_assembly(seq, enzyme, linear=False):
    cuts = re_digest(seq, enzyme, linear=linear)
    return _digest_for_assembly(seq, cuts)


def join_seqs(seqs):
    # sequences = []
    # SeqRecord()
    # every element of seqs could be a Seq or SeqRecord
    # join all annotations
    # join all letter_annotations (intersection of all)
    # assembly = Seq.SeqRecord("", alphabet)
    # assembly = deepcopy(seqs[0][0])
    return seqs


def anneal_oligos():
    # align, find overhangs
    # add feature to seqrecord with name of part (?)
    pass
    # return (overhang1, SeqRecord, overhang2)


def _check_seq_compatibility(seq1, seq2):
    _, overhang1_1, overhang1_2 = seq1
    _, overhang2_1, overhang2_2 = seq2
    return (overhang1_2[0], overhang1_2[1]) == overhang2_1


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


def assemble_sequences(seqs, linear=True):
    alphabet = seqs[0][0].seq.alphabet
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
    # add circular annotation?
    return joined_seq
