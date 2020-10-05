from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition


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


def get_seq(seq):
    if hasattr(seq, "seq"):
        return seq.seq
    else:
        return seq


def join_seqs(seqs):
    to_concat = []
    features = []
    offset = 0
    for seq in seqs:
        to_concat.append(get_seq(seq))
        if hasattr(seq, "features"):
            for feature in seq.features:
                features.append(feature._shift(offset))
        offset += len(seq)
    concatenated_seq = sum(to_concat, Seq(""))
    return SeqRecord(concatenated_seq, features=features)


def anneal_oligos():
    # align, find overhangs
    # add feature to seqrecord with name of part (?)
    pass
