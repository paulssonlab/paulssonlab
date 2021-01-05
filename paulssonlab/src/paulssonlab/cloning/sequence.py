from numbers import Integral
from math import copysign
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition
from paulssonlab.util import sign, format_sign

MAX_SEQUENCE_STR_LENGTH = 34
SEQUENCE_OVERHANG_STR_BUFFER = 4


def _ligate_goldengate(seq1, seq2):
    pass


def _ligate_gibson(seq1, seq2):
    pass


LIGATION_METHODS = {"goldengate": _ligate_goldengate, "gibson": _ligate_gibson}


class DsSeqRecord(SeqRecord):
    def __init__(
        self,
        seq,
        circular=False,
        upstream_overhang=0,
        downstream_overhang=0,
        id="<unknown id>",
        name="<unknown name>",
        description="<unknown description>",
        features=None,
        annotations=None,
        letter_annotations=None,
    ):
        self.circular = circular
        self.upstream_overhang = upstream_overhang or 0
        self.downstream_overhang = downstream_overhang or 0
        if isinstance(seq, SeqRecord):
            id = seq.id or id
            name = seq.name or name
            description = seq.description or description
            features = seq.features.copy() or features
            annotations = seq.annotations.copy() or annotations
            letter_annotations = seq.letter_annotations.copy() or letter_annotations
            seq = seq.seq
        super().__init__(
            seq,
            id=id,
            name=name,
            description=description,
            features=features,
            annotations=annotations,
            letter_annotations=letter_annotations,
        )

    def ligate(self, seq, method="goldengate", try_reverse_complement=True, **kwargs):
        # TAKE seq or seqs as input, join with best and return rest of pool?
        score1, ligated1 = self._ligate(seq, method=method, **kwargs)
        if try_reverse_complement:
            score2, ligated2 = self._ligate(
                reverse_complement(seq), method=method, **kwargs
            )
            # negative score means sequences cannot be ligated (e.g., incompatible sticky stops)
            if max(score1, score2) < 0:
                raise ValueError("attempting to ligate incompatible sequences")
            if score1 >= score2:
                return ligated1
            else:
                return ligated2
        # TODO: raise error if cannot ligate

    def _ligate(self, seq, method="goldengate", **kwargs):
        pass
        # return score

    def annotate(self, label, overhangs=True):
        pass

    def reverse_complement(self):
        rc = self.__class__(super().reverse_complement())
        rc.upstream_overhang, rc.downstream_overhang = (
            -self.downstream_overhang,
            -self.upstream_overhang,
        )
        return rc

    def slice(self, start, stop, annotation_start=None, annotation_stop=None):
        if start is None:
            start = 0
        if stop is None:
            stop = len(seq)
        if start < 0:
            start = start % len(seq)
        if stop < 0:
            stop = stop % len(seq)
        if stop < start:
            if not self.circular:
                raise ValueError(
                    "attempting to slice a non-circular sequence with stop < start"
                )
            slice1 = seq.slice(start, None, annotation_start=annotation_start)
            slice2 = seq.slice(None, stop, annotation_stop=annotation_stop)
            new_seq = slice1 + slice2
            # TODO: join features that wrap around? using a join_features func?
        else:
            new_seq = SeqRecord.__getitem__(seq, slice(start, stop))
            # trim overhangs if necessary
            # new_seq.upstream_overhang = math.copysign(min(abs(start)))
            new_seq.upstream_overhang = copysign(
                max(abs(self.upstream_overhang) - start, 0), self.upstream_overhang
            )
            new_seq.downstream_overhang = copysign(
                max(len(seq) - abs(self.downstream_overhang) - stop, 0),
                self.downstream_overhang,
            )
            # copy features
            features = []
            for feature in seq.features:
                new_feature = _slice_seqrecord_feature(
                    feature,
                    start,
                    stop,
                    annotation_start=annotation_start,
                    annotation_stop=annotation_stop,
                )
                if new_feature is not None:
                    features.append(new_feature)
            new_seq.features = features
            # copy letter_annotations
            # NOTE: annotation_start and annotation_stop only affect features
            # letter_annotations are dropped outside of start:stop
            for key, value in seq.letter_annotations.items():
                new_seq._per_letter_annotations[key] = value[start:stop]
        return new_seq

    def __getitem__(self, index):
        # include truncated annotations? (unlike SeqRecord)
        if isinstance(index, Integral):
            return self.seq[index]
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(self))
            if step != 1:
                raise ValueError("step != 1 not allowed")
            return self.slice(start, stop)
        else:
            raise ValueError("invalid index")

    def __add__(self, other):
        self_overhang_seq = self.seq[-abs(self.downstream_overhang) :]
        other_overhang_seq = other.seq[: abs(other.upstream_overhang)]
        if not (
            self.downstream_overhang == -other.upstream_overhang
            and self_overhang_seq.lower() == other_overhang_seq.lower()
        ):
            raise ValueError(
                f"attempting to join incompatible overhangs: {self.downstream_overhang_seq} ({format_sign(self.downstream_overhang)}) with {other.upstream_overhang_seq} ({format_sign(other.upstream_overhang)})"
            )
        # trim overhang
        other = other[abs(self.downstream_overhang) :]
        new = self.__class__(
            self.seq + other.seq,
            upstream_overhang=self.upstream_overhang,
            downstream_overhang=other.downstream_overhang,
        )
        length = len(self)
        for f in other.features:
            new.features.append(f._shift(length))
        if self.id == other.id:
            new.id = self.id
        if self.name == other.name:
            new.name = self.name
        if self.description == other.description:
            new.description = self.description
        for k, v in self.annotations.items():
            if k in other.annotations and other.annotations[k] == v:
                new.annotations[k] = v
        for k, v in self.letter_annotations.items():
            if k in other.letter_annotations:
                new.letter_annotations[k] = v + other.letter_annotations[k]
        return new

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(seq={self.seq!r},"
            f" upstream_overhang={self.upstream_overhang}, downstream_overhang={self.downstream_overhang},"
            f" id={self.id!r},"
            f" name={self.name!r}, description={self.description!r})"
        )

    def __str__(self):
        seq = str(self.seq)
        if hasattr(self.seq, "complement"):
            seq_rc = str(self.seq.complement())
        else:
            seq_rc = seq
        if len(seq) > 80:
            upstream_length = max(
                abs(self.upstream_overhang) + SEQUENCE_OVERHANG_STR_BUFFER,
                MAX_SEQUENCE_STR_LENGTH,
            )
            downstream_length = max(
                abs(self.downstream_overhang) + SEQUENCE_OVERHANG_STR_BUFFER,
                MAX_SEQUENCE_STR_LENGTH,
            )
            seq = seq[:upstream_length] + "..." + seq[-downstream_length:]
            seq_rc = seq_rc[:upstream_length] + "..." + seq_rc[-downstream_length:]
        lines = []
        for strand in (-1, 1):
            if strand == 1:
                formatted_seq = seq_rc
            else:
                formatted_seq = seq
            if sign(self.upstream_overhang) == strand:
                length = abs(self.upstream_overhang)
                formatted_seq = " " * length + formatted_seq[length:]
            if sign(self.downstream_overhang) == strand:
                length = abs(self.downstream_overhang)
                formatted_seq = formatted_seq[:-length] + " " * length
            lines.append(formatted_seq)
        # TODO: do we want to output this metadata?
        # if self.id:
        #     lines.append(f"ID: {self.id}")
        # if self.name:
        #     lines.append(f"Name: {self.name}")
        # if self.description:
        #     lines.append(f"Description: {self.description}")
        # lines.append(f"Number of features: {len(self.features)}")
        # for a in self.annotations:
        #     lines.append(f"/{a}={str(self.annotations[a])}")
        return "\n".join(lines)


def ligate(seqs, method="goldengate"):
    # handle seq/seqrecord/dsseqrecord seqs by promoting all to dsseqrecords
    pass


def join_seqs(seqs):
    to_concat = []
    features = []
    offset = 0
    for seq in seqs:
        to_concat.appstop(get_seq(seq))
        if hasattr(seq, "features"):
            for feature in seq.features:
                features.appstop(feature._shift(offset))
        offset += len(seq)
    concatenated_seq = sum(to_concat, Seq(""))
    return SeqRecord(concatenated_seq, features=features)


def reverse_complement(seq):
    if hasattr(seq, "features"):
        return seq.reverse_complement(
            id=True, name=True, description=True, annotations=True, dbxrefs=True
        )
    elif hasattr(seq, "reverse_complement"):
        return seq.reverse_complement()
    else:
        return Seq(seq).reverse_complement()


def _slice_seqrecord_feature(
    feature, start, stop, annotation_start=None, annotation_stop=None
):
    if annotation_start is None:
        annotation_start = start
    if annotation_stop is None:
        annotation_stop = stop
    if (
        annotation_stop <= feature.location.nofuzzy_start
        or annotation_start >= feature.location.nofuzzy_end
    ):
        return None
    new_feature = feature._shift(-start)
    start_loc = new_feature.location.start
    stop_loc = new_feature.location.end
    if annotation_start > feature.location.nofuzzy_start:
        start_loc = ExactPosition(start - annotation_start)
    if annotation_stop < feature.location.nofuzzy_end:
        stop_loc = ExactPosition(annotation_stop - start)
    new_feature.location = FeatureLocation(
        start_loc, stop_loc, strand=new_feature.location.strand
    )
    return new_feature


def slice_seq(seq, start, stop, annotation_start=None, annotation_stop=None):
    if isinstance(seq, DsSeqRecord):
        return seq.slice(
            start,
            stop,
            annotation_start=annotation_start,
            annotation_stop=annotation_stop,
        )
    if start is None:
        start = 0
    if stop is None:
        stop = len(seq)
    if start < 0:
        start = start % len(seq)
    if stop < 0:
        stop = stop % len(seq)
    if stop < start:
        slice1 = slice_seq(seq, start, None, annotation_start=annotation_start)
        slice2 = slice_seq(seq, None, stop, annotation_stop=annotation_stop)
        new_seq = slice1 + slice2
        # TODO: join features that wrap around? using a join_features func?
    else:
        new_seq = seq[start:stop]
        if hasattr(seq, "features"):
            features = []
            for feature in seq.features:
                new_feature = _slice_seqrecord_feature(
                    feature,
                    start,
                    stop,
                    annotation_start=annotation_start,
                    annotation_stop=annotation_stop,
                )
                if new_feature is not None:
                    features.append(new_feature)
            new_seq.features = features
        if hasattr(seq, "letter_annotations"):
            # annotation_start and annotation_stop only affect features
            # letter_annotations are dropped outside of start:stop
            for key, value in seq.letter_annotations.items():
                new_seq._per_letter_annotations[key] = value[start:stop]
    return new_seq


def get_seq(seq):
    if hasattr(seq, "seq"):
        return seq.seq
    else:
        return seq


def smoosh_sequences(a, b):
    """Returns the shortest sequence beginning with `a` and stoping with `b`.

    For example, this will eliminate a common substring on the junction
    between the two sequences.
    """
    min_length = min(len(a), len(b))
    for i in range(min_length):
        a_tail = a[-min_length + i :]
        b_head = b[: min_length - i]
        if a_tail == b_head:
            return a + b[min_length - i :]
    return a + b
