import re
from dataclasses import dataclass
from enum import Enum

import Bio.Restriction


class CutDirection(Enum):
    UPSTREAM = 1
    DOWNSTREAM = 2
    BOTH = 3

    def __add__(self, other):
        if (self == self.UPSTREAM and other == self.DOWNSTREAM) or (
            self == self.DOWNSTREAM and other == self.UPSTREAM
        ):
            return self.BOTH
        else:
            raise ValueError(f"cannot combine {self} and {other}")

    def reverse(self):
        if self == self.BOTH:
            return self
        elif self == self.UPSTREAM:
            return self.DOWNSTREAM
        elif self == self.DOWNSTREAM:
            return self.UPSTREAM
        else:
            raise ValueError(f"cannot reverse {self}")


@dataclass
class CutSite:
    cut5: int
    cut3: int
    cut_direction: CutDirection


def _re_search(seq, enzyme, circular=None):
    # TODO: hack to avoid circular deps
    from paulssonlab.cloning.sequence import get_seq

    seq = str(get_seq(seq))
    compsite = re.compile(
        enzyme.compsite.pattern, enzyme.compsite.flags | re.IGNORECASE
    )
    if circular:
        seq = seq + seq[: enzyme.size - 1]
    re_sites = [
        (
            i.start(),
            CutDirection.UPSTREAM
            if i.group(1) is not None
            else CutDirection.DOWNSTREAM,
        )
        for i in re.finditer(compsite, seq)
    ]
    return re_sites


def _re_search_cuts(binding_locs, enzyme, length):
    cuts = {}
    for loc, cut_direction in binding_locs:
        for cut5, cut3 in ((enzyme.fst5, enzyme.fst3), (enzyme.scd5, enzyme.scd3)):
            if cut5 is None and cut3 is None:
                continue
            if cut_direction == CutDirection.UPSTREAM:
                if cut5 is not None:
                    cut5_loc = loc + cut5
                else:
                    cut5_loc = None
                if cut3 is not None:
                    cut3_loc = loc + enzyme.size + cut3
                else:
                    cut3_loc = None
            elif cut_direction == CutDirection.DOWNSTREAM:
                if cut3 is not None:
                    cut5_loc = loc - cut3
                else:
                    cut5_loc = None
                if cut5 is not None:
                    cut3_loc = loc + enzyme.size - cut5
                else:
                    cut3_loc = None
            else:
                raise ValueError(f"not expecting cut direction {cut_direction}")
            cut = CutSite(cut5_loc, cut3_loc, cut_direction)
            key = (cut5_loc % length, cut3_loc % length)
            # uniquify
            if key in cuts:
                old_cut = cuts[key]
                old_cut.cut_direction += cut_direction
            else:
                cuts[key] = cut
    return list(cuts.values())


def re_search(seq, enzyme, circular=None):
    if hasattr(seq, "circular") and circular is None:
        circular = seq.circular
    if isinstance(enzyme, str):
        enzyme = getattr(Bio.Restriction, enzyme)
    binding_locs = _re_search(seq, enzyme, circular=circular)
    cuts = _re_search_cuts(binding_locs, enzyme, len(seq))
    return cuts


def _re_digest(seq, cuts):
    offset = 0
    seqs = []
    for cut in cuts:
        frags, new_offset = seq.cut(cut.cut5 - offset, cut.cut3 - offset)
        if len(frags) == 1:
            frags[0].upstream_cut_direction = cut.cut_direction
            frags[0].downstream_cut_direction = cut.cut_direction.reverse()
        elif len(frags) == 2:
            frags[0].downstream_cut_direction = cut.cut_direction.reverse()
            frags[1].upstream_cut_direction = cut.cut_direction
        else:
            raise NotImplementedError
        seqs.extend(frags[:-1])
        seq = frags[-1]
        offset += new_offset
    seqs.append(seq)
    return seqs


def re_digest(seq, enzyme, circular=None):
    # TODO: hack to avoid circular deps
    from paulssonlab.cloning.sequence import DsSeqRecord

    if not isinstance(seq, DsSeqRecord):
        seq = DsSeqRecord(seq, circular=circular)
    cuts = re_search(seq, enzyme, circular=circular)
    return _re_digest(seq, cuts)
