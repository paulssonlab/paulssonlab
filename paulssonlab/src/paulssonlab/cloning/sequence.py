from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from math import ceil
from numbers import Integral
from cytoolz import partial
from itertools import product as it_product
from collections import defaultdict, OrderedDict
from operator import attrgetter
from copy import deepcopy
from typing import Union, NamedTuple
from paulssonlab.cloning.enzyme import re_digest
from paulssonlab.util import sign, format_sign


MAX_SEQUENCE_STR_LENGTH = 34
SEQUENCE_OVERHANG_STR_BUFFER = 4


def get_seq(seq):
    if hasattr(seq, "seq"):
        return seq.seq
    else:
        return seq


def reverse_complement(seq):
    if seq is None:
        raise ValueError("cannot reverse complement None")
    if hasattr(seq, "features"):
        return seq.reverse_complement(
            id=True, name=True, description=True, annotations=True, dbxrefs=True
        )
    elif hasattr(seq, "reverse_complement"):
        return seq.reverse_complement()
    else:
        return Seq(seq).reverse_complement()


# we need this alias for enumerate_matches
reverse_complement_ = reverse_complement


# def replace(seq, old, new):
#     seq = str(get_seq(seq))  # TODO: handle native seq objects (to preserve features)
#     return seq.replace(old, new).replace(
#         sequence.reverse_complement(old), sequence.reverse_complement(new)
#     )


def _assemble_goldengate(
    seq1, seq2, junction_annotation="GG overhang", keep_junction_annotations=True
):
    if seq2 is None:
        seq2 = seq1
        circularizing = True
    else:
        circularizing = False
    if seq1.can_ligate(seq2):
        junction_length = abs(seq1.downstream_overhang)
        if circularizing:
            if junction_length == 0:
                return None, -1
            product = seq1.circularize()
        elif keep_junction_annotations:
            product = seq1 + seq2
        else:
            product = seq1.slice(
                None, None, annotation_stop=len(seq1) - junction_length
            ) + seq2.slice(None, None, annotation_start=junction_length)
        return product, junction_length
    else:
        return None, -1


def _assemble_gibson(
    seq1,
    seq2,
    min_overlap=15,
    max_overlap=200,
    junction_annotation="Gibson homology",
    junction_annotation_type="misc_feature",
    keep_junction_annotations=True,
):
    if seq2 is None:
        seq2 = seq1
        circularizing = True
        max_overlap = min(
            len(seq1) - 1, max_overlap
        )  # prevent junction from being the entire sequence
    else:
        circularizing = False
    # TODO: we ignore overhangs and assume the full sequence is double-stranded
    # this is probably not too unreasonable for gibson
    junction_length = find_homologous_ends(
        seq1.seq_lower(), seq2.seq_lower(), max_overlap=max_overlap
    )
    if junction_length == max_overlap:
        raise ValueError(
            "found Gibson junction with max_overlap, the results may be nonsensical"
        )
    if junction_length is not None:
        if circularizing:
            product = (
                seq1.fill_in()
                .slice(0, len(seq1) - junction_length, annotation_stop=len(seq1))
                .circularize()
            )
        elif keep_junction_annotations:
            product = seq1.fill_in("downstream") + seq2.fill_in("upstream").slice(
                junction_length, None, annotation_start=0
            )
        else:
            product = seq1.fill_in("downstream").slice(
                None, None, annotation_stop=len(seq1) - junction_length
            ) + seq2.fill_in("upstream").slice(junction_length, None)
        if junction_length < min_overlap:
            # instead of returning junction_length and thresholding in assemble,
            # we threshold in this function so that any thresholding kwargs can be handled
            # in a method-dependent way
            return product, -1
        if junction_annotation is not None:
            feature = SeqFeature(
                FeatureLocation(len(seq1) - junction_length, len(seq1)),
                type=junction_annotation_type,
            )
            feature.qualifiers["label"] = [junction_annotation]
            product.features.append(feature)
        return product, junction_length
    else:
        return None, -1


def ensure_dsseqrecords(*seqs):
    seqs = list(seqs)
    for idx in range(len(seqs)):
        if isinstance(seqs[idx], str):
            seqs[idx] = Seq(seqs[idx])
        if not isinstance(seqs[idx], DsSeqRecord):
            seqs[idx] = DsSeqRecord(seqs[idx])
    return seqs


class DsSeqRecord(SeqRecord):
    assembly_methods = {"goldengate": _assemble_goldengate, "gibson": _assemble_gibson}

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
            features = features or deepcopy(seq.features)
            annotations = annotations or deepcopy(seq.annotations)
            letter_annotations = letter_annotations or deepcopy(seq.letter_annotations)
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
            # need to grab the actual sequence
            seq = seq.seq
        elif isinstance(seq, str):
            seq = Seq(seq)
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
        self.annotations.setdefault("molecule_type", "ds-DNA")  # default to ds-DNA
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

    def __bool__(self):
        return bool(len(self))

    @property
    def circular(self):
        return self.annotations["topology"] == "circular"

    @circular.setter
    def circular(self, value):
        self.annotations["topology"] = "circular" if value else "linear"

    def _check_overhang(self, overhang1, overhang2):
        if len(self) and abs(overhang1) + abs(overhang2) >= len(self):
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
        return self.seq_lower()[: abs(self.upstream_overhang)]

    @property
    def downstream_overhang_seq(self):
        return self.seq_lower()[len(self) - abs(self.downstream_overhang) : len(self)]

    def can_ligate(self, other):
        return (
            self.downstream_overhang == -other.upstream_overhang
            and self.downstream_overhang_seq == other.upstream_overhang_seq
        )

    def can_circularize(self):
        if self.circular:
            return True
        return self.can_ligate(self)

    def circularize(self):
        if self.circular:
            return self
        if not self.can_circularize():
            raise ValueError(
                f"attempting to circularize by ligating incompatible overhangs: {self.downstream_overhang_seq} ({format_sign(self.downstream_overhang)}) with {self.upstream_overhang_seq} ({format_sign(self.upstream_overhang)})"
            )
        seq = self[: len(self) - abs(self.downstream_overhang)]
        seq.upstream_overhang = 0
        seq.upstream_inward_cut = None
        seq.circular = True
        return seq

    def assemble(
        self,
        seq=None,
        method="goldengate",
        try_reverse_complement=True,
        circularize=True,
        **kwargs,
    ):
        if seq is not None:
            if not isinstance(seq, self.__class__):
                seq = self.__class__(seq)
            if self.circular or seq.circular:
                raise ValueError("cannot assemble a circular sequence")
            if try_reverse_complement:
                possible_assemblies = it_product(
                    (self, self.reverse_complement()), (seq, reverse_complement(seq))
                )
            else:
                possible_assemblies = []
            results = [
                l[0]._assemble(l[1], method=method, **kwargs)
                for l in possible_assemblies
            ]
            results = sorted(
                [r for r in results if r[1] > 0], key=lambda x: x[1], reverse=True
            )
            # negative score means sequences cannot be ligated (e.g., incompatible sticky ends)
            if not len(results):
                raise ValueError("attempting to assemble incompatible sequences")
            # try circularizing
            product = results[0][0]
        else:
            product = self
        if circularize:
            # you must explicitly circularize if you want to
            # ligate blunt ends
            circularized, score = product._assemble(None, method=method, **kwargs)
            if score > 0:
                product = circularized
        return product

    def _assemble(self, seq, method="goldengate", **kwargs):
        return self.assembly_methods[method](self, seq, **kwargs)

    def annotate(self, label, type="misc_feature", overhangs=True):
        if overhangs:
            start = 0
            stop = len(self)
        else:
            start = abs(self.upstream_overhang)
            stop = len(self) - abs(self.downstream_overhang)
        feature = SeqFeature(FeatureLocation(start, stop), type=type)
        feature.qualifiers["label"] = [label]
        return self.__class__(self, features=[*self.features, feature])

    def annotate_overhangs(self, type="misc_feature"):
        features = []
        if self.upstream_overhang:
            upstream_feature = SeqFeature(
                FeatureLocation(
                    0, abs(self.upstream_overhang), strand=self.upstream_strand
                ),
                type=type,
            )
            upstream_feature.qualifiers["label"] = ["overhang"]
            features.append(upstream_feature)
        if self.downstream_overhang:
            downstream_feature = SeqFeature(
                FeatureLocation(
                    len(self) - abs(self.downstream_overhang),
                    len(self),
                    strand=self.downstream_strand,
                ),
                type=type,
            )
            downstream_feature.qualifiers["label"] = ["overhang"]
            features.append(downstream_feature)
        return self.__class__(self, features=[*features, *self.features])

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

    def lower(self):
        new = self.__class__(self)
        new.seq = new.seq.lower()
        return new

    def seq_lower(self):
        return str(self.seq).lower()

    def reindex(self, loc):
        if not self.circular:
            raise ValueError("cannot reindex a non-circular sequence")
        new = self[loc:] + self[:loc]
        new.circular = True
        return new

    def fill_in(self, ends="both"):
        new = self.__class__(self)
        if ends in ("both", "upstream"):
            new.upstream_overhang = 0
            new.upstream_inward_cut = None
        if ends in ("both", "downstream"):
            new.downstream_overhang = 0
            new.downstream_inward_cut = None
        return new

    def cut(self, cut5, cut3):
        if not (-len(self) <= cut5 <= len(self) and -len(self) <= cut3 <= len(self)):
            raise IndexError("attempting to cut out of bounds")
        if cut5 < 0:
            cut5 += len(self)
        if cut3 < 0:
            cut3 += len(self)
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
                    <= len(self)
                    - abs(self.downstream_overhang) * (self.downstream_overhang < 0)
                )
                and (
                    abs(self.upstream_overhang) * (self.upstream_overhang < 0)
                    <= cut3
                    <= len(self)
                    - abs(self.downstream_overhang) * (self.downstream_overhang > 0)
                )
            ):
                raise ValueError("attempting to cut a missing strand")
            first = self[:max_loc]
            second = self[min_loc:]
            second.upstream_overhang = -overhang
            first.downstream_overhang = overhang
            seqs = [first, second]
            return seqs, min_loc

    def slice_or_reindex(self, start, stop, **kwargs):
        if start == stop:
            return self.reindex(start)
        else:
            return self.slice(start, stop, **kwargs)

    def slice(self, start, stop, annotation_start=True, annotation_stop=True):
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
            # SeqRecord.__getitem__ does not preserve annotations or features
            new_seq = SeqRecord.__getitem__(self, slice(start, stop))
            # we need to copy here, or changing the annotations
            # (e.g., circular) of the copy will change those of the original
            new_seq.annotations = deepcopy(self.annotations)
            # because circular is stored in annotations, need to set this after setting annotations
            new_seq.circular = False
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
        if isinstance(index, Integral):
            return self.seq[index]
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(self))
            if step != 1:
                raise ValueError("step != 1 not allowed")
            return self.slice(start, stop)
        else:
            raise ValueError(f"invalid index: {index}")

    def __add__(self, other):
        if isinstance(other, Seq):
            other = self.__class__(other)
        if self.circular:
            raise ValueError("cannot append to a circular sequence")
        if other.circular:
            raise ValueError("cannot append a circular sequence")
        if not self.can_ligate(other):
            raise ValueError(
                f"attempting to ligate incompatible overhangs: {self.downstream_overhang_seq} ({format_sign(self.downstream_overhang)}) with {other.upstream_overhang_seq} ({format_sign(other.upstream_overhang)})"
            )
        # trim overhang
        # other[abs(self.downstream_overhang) :] would truncate annotations
        other = other.slice(
            abs(self.downstream_overhang),
            None,
            annotation_start=None,
            annotation_stop=None,
        )
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

    def __radd__(self, other):
        if isinstance(other, Seq):
            other = self.__class__(other)
        return other.__add__(self)

    def __eq__(self, other):
        # TODO
        return self.seq.lower() == other.seq.lower()

    def __ne__(self, other):
        return not self.__eq__(other)

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


def assemble(seqs, **kwargs):
    # we need only cast the first seq to DsSeqRecord because DsSeqRecord.assemble casts its input
    if not len(seqs):
        return DsSeqRecord(Seq(""))
    product = seqs[0]
    if not isinstance(product, DsSeqRecord):
        product = DsSeqRecord(product)
    for seq in seqs[1:]:
        product = product.assemble(seq, **kwargs)
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
    feature, start, stop, annotation_start=True, annotation_stop=True
):
    # annotation_start/stop=True means cut annotations along sequence slice boundaries
    # =False causes unintuitive behavior, so set to None instead to never truncate annotations
    if annotation_start is True:
        annotation_start = start
    elif annotation_start is False:
        annotation_start = None
    if annotation_stop is True:
        annotation_stop = stop
    elif annotation_stop is False:
        annotation_stop = None
    if (
        annotation_stop is not None and annotation_stop <= int(feature.location.start)
    ) or (
        annotation_start is not None and annotation_start >= int(feature.location.end)
    ):
        return None
    new_feature = feature._shift(-start)  # copy feature and shift
    start_loc = new_feature.location.start
    stop_loc = new_feature.location.end
    if annotation_start is not None and annotation_start > int(feature.location.start):
        start_loc = ExactPosition(annotation_start - start)
    if annotation_stop is not None and annotation_stop < int(feature.location.end):
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


def find_homologous_ends(a, b, max_overlap=None):
    _max_overlap = min(len(a), len(b))
    if max_overlap is not None:
        max_overlap = min(_max_overlap, max_overlap)
    else:
        max_overlap = _max_overlap
    for overlap in reversed(range(1, max_overlap + 1)):
        a_tail = a[-overlap:]
        b_head = b[:overlap]
        if a_tail == b_head:
            return overlap
    return 0


def smoosh_sequences(*seqs, max_overlap=None):
    seq = seqs[0]
    if len(seqs) >= 2:
        for seq2 in seqs[1:]:
            seq = _smoosh_sequences(seq, seq2)
    return seq


def _smoosh_sequences(a, b, max_overlap=None):
    """Returns the shortest sequence starting with `a` and ending with `b`.

    For example, this will eliminate a common substring on the junction
    between the two sequences.
    """
    overlap = find_homologous_ends(a, b, max_overlap=max_overlap)
    if overlap is not None:
        return a + b[overlap:]
    else:
        return a + b


def find_aligned_subsequence(seq, subseq, last=False):
    period = len(subseq)
    idxs = range(int(ceil(len(seq) / period)))
    if last:
        idxs = reversed(idxs)
    for i in idxs:
        if seq[i * period : (i + 1) * period] == subseq:
            return i * period
    return None


def count_matching(a, b):
    max_length = min(len(a), len(b))
    score = 0
    for i in range(max_length):
        if a[i] == b[i]:
            score += 1
    return score, 0, max_length


def longest_justified_matching(a, b, right=False):
    max_length = min(len(a), len(b))
    score = 0
    idxs = range(max_length)
    if right:
        idxs = reversed(idxs)
    for i in idxs:
        if a[i] == b[i]:
            score += 1
        else:
            break
    if right:
        start = max_length - score
        stop = max_length
    else:
        start = 0
        stop = score
    return score, start, stop


def longest_contiguous_matching(a, b):
    max_length = min(len(a), len(b))
    max_score = 0
    start = 0
    score = 0
    idxs = range(max_length)
    for i in idxs:
        if a[i] == b[i]:
            score += 1
        else:
            if score > max_score:
                start = i - score
                max_score = score
            score = 0
    if score > max_score:
        start = max_length - score
        max_score = score
    stop = start + max_score
    return max_score, start, stop


def iterate_shifts(a, b):
    if a.circular is False and b.circular is False:
        for shift in range(1, len(a) + len(b)):
            a_start = max(shift - len(b), 0)
            a_stop = min(shift, len(a))
            b_start = a_start + len(b) - shift
            b_stop = a_stop + len(b) - shift
            yield a_start, a_stop, b_start, b_stop
    elif a.circular is False and b.circular is True:
        for b_start, b_stop, a_start, a_stop in iterate_shifts(b, a):
            yield a_start, a_stop, b_start, b_stop
    elif a.circular is True and b.circular is False:
        max_shift = max(len(b) - len(a), 0)
        for shift in range(max_shift + 1):
            # remember: max_shift=shift=0 if len(b) < len(a)
            b_start = max_shift - shift
            b_stop = len(b) - shift
            for a_start in range(len(a)):
                a_stop = (b_stop - b_start + a_start) % len(a)
                yield a_start, a_stop, b_start, b_stop
    elif a.circular is True and b.circular is True:
        raise ValueError("cannot bind a circular sequence to another circular sequence")


class SeqMatch(NamedTuple):
    strand: int
    score: int
    seq1_start: int
    seq1_stop: int
    seq2_start: int
    seq2_stop: int


def iter_matches(
    seq1,
    seq2,
    reverse_complement=None,
    scoring_func=count_matching,
    min_score=10,
    lower=True,
):
    orig_seq1 = seq1
    if reverse_complement is None:
        strands = (1, -1)
    elif reverse_complement is True:
        strands = (-1,)
    elif reverse_complement is False:
        strands = (1,)
    else:
        raise ValueError("expecting reverse_complement to be one of: True, False, None")
    sites = []
    for strand in strands:
        if strand == -1:
            seq1 = reverse_complement_(orig_seq1)
        else:
            seq1 = orig_seq1
        if lower:
            seq1_lower = seq1.lower()
            seq2_lower = seq2.lower()
        else:
            seq1_lower = seq1
            seq2_lower = seq2
        for (seq1_start, seq1_stop, seq2_start, seq2_stop) in iterate_shifts(
            seq1, seq2
        ):
            seq1_overlap = seq1_lower.slice_or_reindex(seq1_start, seq1_stop)
            seq2_overlap = seq2_lower.slice_or_reindex(seq2_start, seq2_stop)
            score, overlap_start, overlap_stop = scoring_func(
                seq1_overlap, seq2_overlap
            )
            if score >= min_score:
                seq1_match_start = seq1_start + overlap_start
                seq2_match_start = seq2_start + overlap_start
                seq1_match_stop = seq1_start + overlap_stop
                seq2_match_stop = seq2_start + overlap_stop
                if strand == -1:
                    # when reverse complementing seq1, need to adjust locations
                    seq1_match_start = len(seq1) - seq1_match_start
                    seq1_match_stop = len(seq1) - seq1_match_stop
                match = SeqMatch(
                    strand,
                    score,
                    seq1_match_start,
                    seq1_match_stop,
                    seq2_match_start,
                    seq2_match_stop,
                )
                yield match


def enumerate_matches(
    seq1,
    seq2,
    reverse_complement=None,
    scoring_func=count_matching,
    min_score=10,
    lower=True,
):
    seq1, seq2 = ensure_dsseqrecords(seq1, seq2)
    sites = iter_matches(
        seq1,
        seq2,
        reverse_complement=reverse_complement,
        scoring_func=scoring_func,
        min_score=min_score,
        lower=lower,
    )
    return sorted(list(sites), key=attrgetter("seq1_start"))


def iter_primer_binding_sites(
    template,
    primer,
    reverse_complement=None,
    scoring_func=None,
    min_score=10,
    require_3prime_clamp=True,
    lower=True,
    return_sequences=False,
):
    template, primer = ensure_dsseqrecords(template, primer)
    if scoring_func is None:
        if require_3prime_clamp is True:
            # require matching at 3' end
            scoring_func = partial(longest_justified_matching, right=True)
        else:
            # do not penalize mismatches
            scoring_func = count_matching
    sites = iter_matches(
        template,
        primer,
        reverse_complement=reverse_complement,
        scoring_func=scoring_func,
        min_score=min_score,
        lower=lower,
    )
    for site in sites:
        if not require_3prime_clamp or site.seq2_stop == len(primer):
            yield site


def enumerate_primer_binding_sites(
    template,
    primer,
    reverse_complement=None,
    scoring_func=None,
    min_score=10,
    require_3prime_clamp=True,
    lower=True,
):
    sites = iter_primer_binding_sites(
        template,
        primer,
        reverse_complement=reverse_complement,
        scoring_func=scoring_func,
        min_score=min_score,
        require_3prime_clamp=require_3prime_clamp,
        lower=lower,
    )
    return sorted(list(sites), key=attrgetter("seq1_start"))


def find_subsequence(template, product, min_score=10):
    sites = iter_matches(
        template, product, scoring_func=longest_contiguous_matching, min_score=min_score
    )
    site = sorted(list(sites), key=attrgetter("score"), reverse=True)[0]
    head = product[: site.seq2_start]
    tail = product[site.seq2_stop :]
    return (site.seq1_start, site.seq1_stop, head, tail)


def _amplicon_location(template, primer1, primer2, min_score=10):
    template, primer1, primer2 = ensure_dsseqrecords(template, primer1, primer2)
    both_sites = []
    sites1 = enumerate_primer_binding_sites(
        template, primer1, reverse_complement=None, min_score=min_score
    )
    if len(sites1) != 1:
        raise ValueError(
            f"expecting a unique primer1 binding site, instead found {len(sites1)}"
        )
    site1 = sites1[0]
    sites2 = enumerate_primer_binding_sites(
        template, primer2, reverse_complement=(sites1[0][0] == 1), min_score=min_score
    )
    if len(sites2) != 1:
        raise ValueError(
            f"expecting a unique primer2 binding site, instead found {len(sites2)}"
        )
    site2 = sites2[0]
    if site1.strand == -1:
        site1, site2 = site2, site1
    if site1.strand != -site2.strand:
        raise ValueError("expecting a forward/reverse primer pair")
    # this is needed to handle negative locs (when working with circular sequences)
    # TODO: not tested!
    loc1 = site1.seq1_stop % len(template)
    loc2 = site2.seq1_stop % len(template)
    return loc1, loc2, site1.score, site2.score


def amplicon_location(template, primer1, primer2, min_score=10):
    template, primer1, primer2 = ensure_dsseqrecords(template, primer1, primer2)
    loc1, loc2, _, _ = _amplicon_location(template, primer1, primer2)
    return loc1, loc2


def pcr(template, primer1, primer2, min_score=10):
    template, primer1, primer2 = ensure_dsseqrecords(template, primer1, primer2)
    loc1, loc2, len1, len2 = _amplicon_location(template, primer1, primer2)
    amplicon = template.slice(
        loc1, loc2, annotation_start=loc1 - len1, annotation_stop=loc2 + len2
    )
    product = primer1 + amplicon + reverse_complement(primer2)
    return product


def anneal(a, b, min_score=6):
    b = reverse_complement(b)
    best_score = -1
    for a_start, a_stop, b_start, b_stop in iterate_shifts(a, b):
        a_overlap = a[a_start:a_stop]
        b_overlap = b[b_start:b_stop]
        score = len(a_overlap)
        if a_overlap == b_overlap and score > best_score:
            best_score = score
            if b_start > 0:
                seq = b[:b_stop] + a[a_stop:]
                upstream_overhang = -b_start
                downstream_overhang = len(a) - a_stop
            else:
                seq = a[:a_stop] + b[b_stop:]
                upstream_overhang = a_start
                downstream_overhang = -(len(b) - b_stop)
    if best_score < min_score:
        raise ValueError(f"could not anneal sequences {a!r} and {b!r}")
    return DsSeqRecord(
        seq,
        upstream_overhang=upstream_overhang,
        downstream_overhang=downstream_overhang,
    )
