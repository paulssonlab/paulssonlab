from typing import NamedTuple, Optional
from paulssonlab.cloning.sequence import reverse_complement, find_primer_binding_site
import paulssonlab.cloning.thermodynamics as thermodynamics
from paulssonlab.cloning.viennarna import dna_secondary_structure, dna_heterodimer


class Primer(NamedTuple):
    seq: str
    binding: str
    tail: str
    tm: Optional[float]
    mfe_monomer: Optional[float]
    mfe_homodimer: Optional[float]

    # @property
    # def length(self):
    #     return len(self.seq)


class PrimerPair(NamedTuple):
    primer1: Primer
    primer2: Primer
    ta: Optional[float]
    mfe_heterodimer: Optional[float]

    @property
    def min_length(self):
        return min(len(self.primer1.seq), len(self.primer2.seq))

    @property
    def max_length(self):
        return max(len(self.primer1.seq), len(self.primer2.seq))


def evaluate_primer(seq, tm_kwargs=None, ta_kwargs=None):
    tm = thermodynamics.tm(seq, **(tm_kwargs or {}))
    ta = thermodynamics.ta_from_tms(tm, **(ta_kwargs or {}))
    mfe_monomer, mfe_homodimer = dna_secondary_structure(seq)
    return Primer(
        seq=seq, tm=tm, ta=ta, mfe_monomer=mfe_monomer, mfe_homodimer=mfe_homodimer
    )


def evaluate_primer_pair(primer1, primer2, tm_kwargs=None, ta_kwargs=None):
    if not isinstance(primer1, Primer):
        primer1 = evaluate_primer(primer1, tm_kwargs=tm_kwargs, ta_kwargs=ta_kwargs)
    if not isinstance(primer2, Primer):
        primer2 = evaluate_primer(primer2, tm_kwargs=tm_kwargs, ta_kwargs=ta_kwargs)
    ta = thermodynamics.ta_from_tms(primer1.tm, primer2.tm, **(ta_kwargs or {}))
    mfe_heterodimer = dna_heterodimer(primer1.seq, primer2.seq)
    return PrimerPair(primer1=primer1, primer2=primer2, mfe_heterodimer=mfe_heterodimer)


def enumerate_primers(
    seq,
    tail="",
    anchor_3prime=False,
    min_length=10,
    max_length=None,
    min_tm=60,
    max_tm=72,
    min_mfe=-8,
    tm_func=thermodynamics.tm,
    monotonic_mfe=False,
    tm_kwargs=None,
):
    primers = []
    assert len(seq) >= 1
    if anchor_3prime:
        stop_idxs = [len(seq) - 1]
    else:
        stop_idxs = reversed(range(len(seq) + 1))
    for stop in stop_idxs:
        if max_length is None:
            min_start = 0
        else:
            min_start = stop - max_length
        for start in reversed(range(min_start, stop - min_length + 1)):
            binding_seq = seq[start:stop]
            primer_seq = tail + binding_seq
            tm = tm_func(binding_seq, **(tm_kwargs or {}))
            if tm < min_tm:
                continue
            elif tm > max_tm:
                break
            mfe_monomer, mfe_homodimer = viennarna.dna_secondary_structure(primer_seq)
            if min(mfe_monomer, mfe_homodimer) < min_mfe:
                if monotonic_mfe:
                    # if we assume mfe monotonically decreases with length, then we can stop extending
                    # this is probably not exactly true but pretty close, and likely speeds things up
                    break
                else:
                    continue
            primer = Primer(
                seq=primer_seq,
                length=len(primer_seq),
                tm=tm,
                mfe_monomer=mfe_monomer,
                mfe_homodimer=mfe_homodimer,
            )
            primers.append(primer)
    return primers


def replace_primer_tail(template, primer, new_tail=""):
    sites = find_primer_binding_site(template, primer, return_sequences=True)
    print(">", sites)
    if len(sites) != 1:
        raise ValueError(
            f"expecting a unique primer binding site, instead got {len(sites)}"
        )
    sense = sites[0][0]  # TODO: NamedTuple-ize
    score = sites[0][2]
    primer_overlap = sites[0][4]
    primer_binding = primer_overlap[len(primer_overlap) - score :]
    if sense == -1:
        new_tail = reverse_complement(new_tail)
    new_primer = new_tail + primer_binding
    return new_primer
