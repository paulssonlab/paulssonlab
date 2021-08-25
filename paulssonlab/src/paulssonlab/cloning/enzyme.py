import re
from typing import NamedTuple
import Bio.Restriction


class CutSite(NamedTuple):
    cut5: int
    cut3: int
    cut_upstream: bool
    recognition_upstream: int
    recognition_downstream: int


def _re_search(seq, enzyme, circular=None):
    # TODO: hack to avoid circular deps
    from paulssonlab.cloning.sequence import get_seq

    seq = str(get_seq(seq))
    compsite = re.compile(
        enzyme.compsite.pattern, enzyme.compsite.flags | re.IGNORECASE
    )
    if circular:
        seq = seq + seq[: enzyme.size - 1]
    re_sites = [(i.start(), i.group(1) is not None) for i in re.finditer(compsite, seq)]
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
            binding_upstream = 0
            binding_downstream = 0
            cuts.append(
                CutSite(
                    cut5_loc,
                    cut3_loc,
                    cut_upstream,
                    binding_upstream,
                    binding_downstream,
                )
            )
    return cuts


def re_search(seq, enzyme, circular=None):
    if hasattr(seq, "circular") and circular is None:
        circular = seq.circular
    if isinstance(enzyme, str):
        enzyme = getattr(Bio.Restriction, enzyme)
    binding_locs = _re_search(seq, enzyme, circular=circular)
    cuts = _re_search_cuts(binding_locs, enzyme)
    return cuts


def _re_digest(seq, cuts):
    offset = 0
    seqs = []
    for cut in cuts:
        frags, new_offset = seq.cut(cut.cut5 - offset, cut.cut3 - offset)
        if len(frags) == 1:
            frags[0].upstream_inward_cut = cut.cut_upstream
            frags[0].downstream_inward_cut = not cut.cut_upstream
        elif len(frags) == 2:
            frags[0].downstream_inward_cut = not cut.cut_upstream
            frags[1].upstream_inward_cut = cut.cut_upstream
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
