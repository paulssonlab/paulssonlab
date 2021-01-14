from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from numbers import Integral
from cytoolz import partial
from itertools import product
from collections import defaultdict, OrderedDict
from paulssonlab.cloning.enzyme import re_digest
from paulssonlab.util import sign, format_sign

MAX_SEQUENCE_STR_LENGTH = 34
SEQUENCE_OVERHANG_STR_BUFFER = 4


def _ligate_goldengate(seq1, seq2):
    pass


def _ligate_gibson(seq1, seq2):
    pass


class DsSeqRecord(SeqRecord):
    ligation_methods = {"goldengate": _ligate_goldengate, "gibson": _ligate_gibson}

    def __init__(
        self,
        seq,
        circular=None,
        upstream_overhang=None,
        downstream_overhang=None,
        upstream_inward_cut=None,
        downstream_inward_cut=None,
        id=None,
        name=None,
        description=None,
        features=None,
        annotations=None,
        letter_annotations=None,
    ):
        self._upstream_overhang = 0
        self._downstream_overhang = 0
        if isinstance(seq, (DsSeqRecord, SeqRecord)):
            id = id or seq.id
            name = name or seq.name
            description = description or seq.description
            features = features or seq.features.copy()
            annotations = annotations or seq.annotations.copy()
            letter_annotations = letter_annotations or seq.letter_annotations.copy()
            seq = seq.seq
            if isinstance(seq, DsSeqRecord):
                circular = circular if circular is not None else seq.circular
                upstream_overhang = (
                    upstream_overhang
                    if upstream_overhang is not None
                    else seq.upstream_overhang
                )
                downstream_overhang = (
                    downstream_overhang
                    if downstream_overhang is not None
                    else seq.downstream_overhang
                )
                upstream_inward_cut = (
                    upstream_inward_cut
                    if upstream_inward_cut is not None
                    else seq.upstream_inward_cut
                )
                downstream_inward_cut = (
                    downstream_inward_cut
                    if downstream_inward_cut is not None
                    else seq.downstream_inward_cut
                )
        id = id or "<unknown id>"
        name = name or "<unknown name>"
        description = description or "<unknown description>"
        super().__init__(
            seq,
            id=id,
            name=name,
            description=description,
            features=features,
            annotations=annotations,
            letter_annotations=letter_annotations,
        )
        if circular is not None:
            self.circular = circular
        elif "topology" not in self.annotations:
            self.circular = False  # default to linear
        if self.circular:
            upstream_overhang = 0
            downstream_overhang = 0
        self.upstream_overhang = upstream_overhang or 0
        self.downstream_overhang = downstream_overhang or 0
        self.upstream_inward_cut = upstream_inward_cut
        self.downstream_inward_cut = downstream_inward_cut

    @property
    def circular(self):
        return self.annotations["topology"] == "circular"

    @circular.setter
    def circular(self, value):
        self.annotations["topology"] = "circular" if value else "linear"

    def _check_overhang(self, overhang1, overhang2):
        if sign(overhang1) == sign(overhang2) != 0:
            if abs(overhang1) + abs(overhang2) > len(self):
                raise ValueError("invalid overhang length")

    @property
    def upstream_overhang(self):
        return self._upstream_overhang

    @upstream_overhang.setter
    def upstream_overhang(self, value):
        self._check_overhang(value, self.downstream_overhang)
        self._upstream_overhang = value

    @property
    def downstream_overhang(self):
        return self._downstream_overhang

    @downstream_overhang.setter
    def downstream_overhang(self, value):
        self._check_overhang(self.upstream_overhang, value)
        self._downstream_overhang = value

    @property
    def upstream_strand(self):
        return sign(self.upstream_overhang)

    @property
    def downstream_strand(self):
        return sign(self.downstream_overhang)

    @property
    def upstream_overhang_seq(self):
        return str(self.seq[: abs(self.upstream_overhang)]).lower()

    @property
    def downstream_overhang_seq(self):
        return str(
            self.seq[len(self) - abs(self.downstream_overhang) : len(self)]
        ).lower()

    def can_ligate(self, other):
        return (
            self.downstream_overhang == -other.upstream_overhang
            and self.downstream_overhang_seq == other.upstream_overhang_seq
        )

    def ligate(self, seq, method="goldengate", try_reverse_complement=True, **kwargs):
        if isinstance(seqs, (Seq, SeqRecord, DsSeqRecord)):
            seqs = [seqs]
        else:
            # cast all sequences to DsSeqRecord if needed
            seqs = [
                DsSeqRecord(s) if not isinstance(s, DsSeqRecord) else s for s in seqs
            ]
        if try_reverse_complement:
            possible_ligations = product(
                (self, self.reverse_complement()), (seq, reverse_complement(seq))
            )
        else:
            possible_ligations = []
        results = [l[0]._ligate(l[1], method=method, **kwargs)]
        results = sorted(
            [r for r in results if r[1] > 0], key=lambda x: x[1], reverse=True
        )
        # negative score means sequences cannot be ligated (e.g., incompatible sticky ends)
        if not len(results):
            raise ValueError("attempting to ligate incompatible sequences")
        return results[0][0]

    def _ligate(self, seq, method="goldengate", **kwargs):
        pass
        # return score
        # if cannot ligate
        return None, -1

    def annotate(self, label, type="misc_feature", overhangs=True):
        if overhangs:
            start = 0
            stop = len(self)
        else:
            start = abs(self.upstream_overhang)
            stop = len(self) - abs(self.downstream_overhang)
        feature = SeqFeature(FeatureLocation(start, stop), type=type)
        feature.qualifiers["label"] = [label]
        return self.__class__(self, features=[feature, *self.features])

    def reverse_complement(self, **kwargs):
        kwargs = {
            **dict(
                id=True, name=True, description=True, annotations=True, dbxrefs=True
            ),
            **kwargs,
        }
        return self.__class__(
            super().reverse_complement(**kwargs),
            upstream_overhang=-self.downstream_overhang,
            downstream_overhang=-self.upstream_overhang,
            upstream_inward_cut=self.downstream_inward_cut,
            downstream_inward_cut=self.upstream_inward_cut,
        )

    def reindex(self, loc):
        if not self.circular:
            raise ValueError("cannot reindex a non-circular sequence")
        new = self[loc:] + self[:loc]
        new.circular = True
        return new

    def cut(self, cut5, cut3):
        if not (-len(self) <= cut5 < len(self) and -len(self) <= cut3 < len(self)):
            raise IndexError("attempting to cut out of bounds")
        cut5 = cut5 % len(self)
        cut3 = cut3 % len(self)
        min_loc = min(cut5, cut3)
        max_loc = max(cut5, cut3)
        overhang = cut5 - cut3
        if self.circular:
            wraparound_length = min_loc + len(self) - max_loc
            if wraparound_length < max_loc - min_loc:
                min_loc, max_loc = max_loc, min_loc
                overhang = -sign(cut5 - cut3) * wraparound_length
            # we do self[min_loc:] + self[:min_loc] instead of self.reindex(min_loc)
            # because the latter returns a circular sequence
            new = self[min_loc:] + self[:min_loc] + self[min_loc:max_loc]
            new.upstream_overhang = -overhang
            new.downstream_overhang = overhang
            return [new], min_loc
        else:
            if not (
                (
                    abs(self.upstream_overhang) * (self.upstream_overhang > 0)
                    <= cut5
                    < len(self)
                    - abs(self.downstream_overhang) * (self.downstream_overhang > 0)
                )
                and (
                    abs(self.upstream_overhang) * (self.upstream_overhang < 0)
                    <= cut3
                    < len(self)
                    - abs(self.downstream_overhang) * (self.downstream_overhang < 0)
                )
            ):
                raise ValueError("attempting to cut a missing strand")
            min_loc = min(cut5, cut3)
            max_loc = max(cut5, cut3)
            overhang = cut5 - cut3
            first = self[:max_loc]
            second = self[min_loc:]
            second.upstream_overhang = -overhang
            first.downstream_overhang = overhang
            return [first, second], min_loc

    def slice(self, start, stop, annotation_start=None, annotation_stop=None):
        if start is None:
            start = 0
        if stop is None:
            stop = len(self)
        if start < 0:
            start = start % len(self)
        if stop < 0:
            stop = stop % len(self)
        if stop < start:
            if not self.circular:
                raise ValueError(
                    "attempting to slice a non-circular sequence with stop < start"
                )
            slice1 = self.slice(start, None, annotation_start=annotation_start)
            slice2 = self.slice(None, stop, annotation_stop=annotation_stop)
            new_seq = slice1 + slice2
        else:
            new_seq = SeqRecord.__getitem__(self, slice(start, stop))
            new_seq.upstream_overhang = max(
                abs(self.upstream_overhang) - start, 0
            ) * sign(self.upstream_overhang)
            new_seq.downstream_overhang = max(
                abs(self.downstream_overhang) - (len(self) - stop), 0
            ) * sign(self.downstream_overhang)
            if start == 0:
                new_seq.upstream_inward_cut = self.upstream_inward_cut
            if stop == len(self):
                new_seq.downstream_inward_cut = self.downstream_inward_cut
            # copy features
            features = []
            for feature in self.features:
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
            for key, value in self.letter_annotations.items():
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
        if self.circular:
            raise ValueError("cannot append to a circular sequence")
        if other.circular:
            raise ValueError("cannot append a circular sequence")
        if not self.can_ligate(other):
            raise ValueError(
                f"attempting to join incompatible overhangs: {self.downstream_overhang_seq} ({format_sign(self.downstream_overhang)}) with {other.upstream_overhang_seq} ({format_sign(other.upstream_overhang)})"
            )
        # trim overhang
        other = other[abs(self.downstream_overhang) :]
        new = self.__class__(
            self.seq + other.seq,
            upstream_overhang=self.upstream_overhang,
            downstream_overhang=other.downstream_overhang,
            upstream_inward_cut=self.upstream_inward_cut,
            downstream_inward_cut=other.downstream_inward_cut,
        )
        length = len(self)
        new.features = _join_features(self.features, other.features, length)
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
            f" upstream_inward_cut={self.upstream_inward_cut}, downstream_inward_cut={self.downstream_inward_cut},"
            f" circular={self.circular!r}, id={self.id!r},"
            f" name={self.name!r}, description={self.description!r})"
        )

    def __str__(self):
        seq = str(self.seq)
        if hasattr(self.seq, "complement"):
            seq_compl = str(self.seq.complement())
        else:
            seq_compl = seq
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
            seq_compl = (
                seq_compl[:upstream_length] + "..." + seq_compl[-downstream_length:]
            )
        lines = []
        for strand in (-1, 1):
            if strand == 1:
                formatted_seq = seq_compl
            else:
                formatted_seq = seq
            if sign(self.upstream_overhang) == strand:
                length = abs(self.upstream_overhang)
                formatted_seq = " " * length + formatted_seq[length:]
            if sign(self.downstream_overhang) == strand:
                length = abs(self.downstream_overhang)
                formatted_seq = formatted_seq[:-length] + " " * length
            if self.circular:
                formatted_seq += " >"
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


def ligate(seqs, **kwargs):
    # cast all sequences to DsSeqRecord if needed
    seqs = [DsSeqRecord(s) if not isinstance(s, DsSeqRecord) else s for s in seqs]
    if not len(seqs):
        return DsSeqRecord(Seq(""))
    product = seqs[0]
    for seq in seqs[1:]:
        product = product.ligate(seq, **kwargs)
    return product


def join_seqs(seqs):
    if any(isinstance(seq, DsSeqRecord) for seq in seqs):
        cls = DsSeqRecord
    elif any(hasattr(seq, "feature") for seq in seqs):
        cls = SeqRecord
    else:
        cls = Seq
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
    seq = cls(concatenated_seq)
    if features:
        seq.features = features
    return seq


def reverse_complement(seq):
    if hasattr(seq, "features"):
        return seq.reverse_complement(
            id=True, name=True, description=True, annotations=True, dbxrefs=True
        )
    elif hasattr(seq, "reverse_complement"):
        return seq.reverse_complement()
    else:
        return Seq(seq).reverse_complement()


def _join_features(a, b, length):
    features = []
    to_join = defaultdict(list)
    for feature in a:
        if int(feature.location.end) < length:
            features.append(feature)
        else:
            key = (feature.type, feature.location_operator, feature.id)
            to_join[key].append(feature)
    for other_feature in b:
        if int(other_feature.location.start) > 0:
            features.append(other_feature._shift(length))
        else:
            key = (
                other_feature.type,
                other_feature.location_operator,
                other_feature.id,
            )
            match = False
            if key in to_join:
                for idx, feature in enumerate(to_join[key]):
                    if feature.qualifiers == other_feature.qualifiers:
                        other_feature_shifted = other_feature._shift(length)
                        location = FeatureLocation(
                            feature.location.start,
                            other_feature_shifted.location.end,
                            strand=feature.location.strand,
                        )
                        new_feature = SeqFeature(
                            location=location,
                            type=feature.type,
                            location_operator=feature.location_operator,
                            id=feature.id,
                            qualifiers=OrderedDict(feature.qualifiers.items()),
                        )
                        features.append(new_feature)
                        to_join[key].pop(idx)
                        match = True
                        break
            if not match:
                features.append(other_feature._shift(length))
    for unmatched_features in to_join.values():
        for feature in unmatched_features:
            features.append(feature)
    return features


def _slice_seqrecord_feature(
    feature, start, stop, annotation_start=None, annotation_stop=None
):
    if annotation_start is None:
        annotation_start = start
    if annotation_stop is None:
        annotation_stop = stop
    if annotation_stop <= int(feature.location.start) or annotation_start >= int(
        feature.location.end
    ):
        return None
    new_feature = feature._shift(-start)  # copy feature and shift
    start_loc = new_feature.location.start
    stop_loc = new_feature.location.end
    if annotation_start > int(feature.location.start):
        start_loc = ExactPosition(start - annotation_start)
    if annotation_stop < int(feature.location.end):
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


def count_matching(a, b):
    s = 0
    for i in range(min(len(a), len(b))):
        if a[i] == b[i]:
            s += 1
    return s


def count_contiguous_matching(a, b, right=False):
    s = 0
    idxs = range(min(len(a), len(b)))
    if right:
        idxs = reversed(idxs)
    for i in idxs:
        if a[i] == b[i]:
            s += 1
        else:
            return s
    return s


def iterate_shifts(a, b):
    for shift in range(1, len(a) + len(b)):
        a_start = max(shift - len(b), 0)
        a_stop = min(shift, len(a))
        b_start = a_start + len(b) - shift
        b_stop = a_stop + len(b) - shift
        a_overlap = a[a_start:a_stop]
        b_overlap = b[b_start:b_stop]
        yield shift, a_start, a_stop, b_start, b_stop


def find_primer_binding_site(
    template,
    primer,
    circular=None,
    try_reverse_complement=True,
    scoring_func=partial(count_contiguous_matching, right=True),
    min_score=6,
):
    if hasattr(template, "circular") and circular is None:
        circular = template.circular
    # TODO: circular
    orig_primer = primer
    # orig_template = template
    if try_reverse_complement:
        strands = (1, -1)
    else:
        strands = (1,)
    sites = []
    for strand in strands:
        if strand == -1:
            primer = reverse_complement(orig_primer)
        else:
            primer = orig_primer
        for (
            shift,
            template_start,
            template_stop,
            primer_start,
            primer_stop,
        ) in iterate_shifts(template, primer):
            template_overlap = template[template_start:template_stop]
            primer_overlap = primer[primer_start:primer_stop]
            score = scoring_func(template_overlap, primer_overlap)
            if score >= min_score:
                sites.append((strand, shift, score))
    return sites


def pcr(template, primer1, primer2, circular=None, min_score=6):
    # cast input sequences to DsSeqRecords if needed
    seqs = [template, primer1, primer2]
    for idx in range(len(seqs)):
        if isinstance(seqs[idx], str):
            seqs[idx] = Seq(seqs[idx])
        if not isinstance(seqs[idx], DsSeqRecord):
            seqs[idx] = DsSeqRecord(seqs[idx])
    template, primer1, primer2 = seqs
    both_sites = []
    for primer in (primer1, primer2):
        sites = find_primer_binding_site(
            template,
            primer,
            circular=circular,
            try_reverse_complement=True,
            min_score=min_score,
        )
        if len(sites) != 1:
            raise ValueError(
                f"expecting a unique primer binding site, instead found {len(sites)}"
            )
        both_sites.append(sites[0])
    if both_sites[0][0] == -1:
        both_sites = both_sites[::-1]
    sense1, loc1, len1 = both_sites[0]
    sense2, loc2, len2 = both_sites[1]
    if sense1 != -sense2:
        raise ValueError("expecting a forward/reverse primer pair")
    start = loc1
    stop = len(template) - loc2
    amplicon = template[start:stop]
    product = primer1 + amplicon + primer2
    return product


def anneal(a, b, min_score=6):
    b = reverse_complement(b)
    loc = None
    best_score = -1
    for shift, a_start, a_stop, b_start, b_stop in iterate_shifts(a, b):
        a_overlap = a[a_start:a_stop]
        b_overlap = b[b_start:b_stop]
        score = len(a_overlap)
        if a_overlap == b_overlap and score > best_score:
            loc = len(b) - shift
            best_score = score
            if b_start > 0:
                seq = b[:b_stop] + a[a_stop:]
                upstream_overhang = -b_start
                downstream_overhang = len(a) - a_stop
            else:
                seq = a[:a_stop] + b[b_stop:]
                upstream_overhang = a_start
                downstream_overhang = -(len(b) - b_stop)
    if loc is None or best_score < min_score:
        raise ValueError(f"could not anneal sequences {a!r} and {b!r}")
    return DsSeqRecord(
        seq,
        upstream_overhang=upstream_overhang,
        downstream_overhang=downstream_overhang,
    )
