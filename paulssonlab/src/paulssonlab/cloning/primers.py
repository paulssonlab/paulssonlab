from functools import cached_property
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from paulssonlab.cloning.sequence import reverse_complement, find_primer_binding_site
import paulssonlab.cloning.thermodynamics as thermodynamics
from paulssonlab.cloning.viennarna import dna_secondary_structure, dna_heterodimer
from paulssonlab.util import any_not_none

# tm_func, tm_parameters
# gc


class Primer:
    #     seq: str = None
    #     binding: str = None
    #     overhang: str = None
    #     tm: Optional[float] = None
    #     mfe_monomer: Optional[float] = None
    #     mfe_homodimer: Optional[float] = None
    #     primer3: dict = field(repr=False, default_factory=dict)

    def __init__(
        self,
        overhang=None,
        binding=None,
        template=None,
        binding_length=None,
        overhang_length=None,
        mfe_monomer=None,
        mfe_homodimer=None,
        primer3=None,
        tm_func=thermodynamics.tm,
        tm_kwargs=None,
    ):
        if mfe_monomer is not None:
            self.mfe_monomer = mfe_monomer
        if mfe_homodimer is not None:
            self.mfe_homodimer = mfe_homodimer
        self.tm_func = tm_func
        self.tm_kwargs = tm_kwargs
        if primer3 is not None:
            if any_not_none(
                overhang, binding, template, binding_length, overhang_length
            ):
                raise ValueError(
                    "if primer3 is specified, cannot also specify other primer sequences/lengths"
                )
            self.primer3 = primer3
            self.binding = Seq(primer3["SEQUENCE"])
            self.overhang = Seq(primer3["OVERHANG"])
            self.gc = primer3["GC_PERCENT"]
        else:
            if overhang is not None:
                if binding is not None:
                    if any_not_none(template, binding_length, overhang_length):
                        raise ValueError(
                            "if binding and overhang are specified, cannot also specify template, binding_length, or overhang_length"
                        )
                    self.overhang = overhang
                    self.binding = binding
                else:
                    if binding_length is not None:
                        if overhang_length is not None:
                            raise ValueError(
                                "cannot specify both binding_length and overhang_length"
                            )
                        overhang_length = len(seq) - binding_length
                        self.binding = overhang[overhang_length:]
                        self.overhang = overhang[:overhang_length]
                    else:
                        if overhang_length is not None:
                            self.binding = overhang[overhang_length:]
                            self.overhang = overhang[:overhang_length]
                        else:
                            if template is not None:
                                overhang, binding, _ = get_primer_overhang(
                                    template, overhang
                                )
                            else:
                                # cannot have a primer with overhang but no binding
                                # so if only one sequence is specified as a positional arg,
                                # assume it is the binding region
                                self.binding = overhang
                                self.overhang = overhang[:0]  # zero-length of same type
            else:
                if binding is not None:
                    self.binding = binding
                    self.overhang = binding[:0]  # zero-length of same type
                else:
                    raise ValueError("need to specify primer sequence")

    @property
    def seq(self):
        return self.overhang + self.binding

    @cached_property
    def gc(self):
        return GC(self.seq)

    @cached_property
    def tm(self):
        return self.tm_func(self.binding, **(self.tm_kwargs or {}))

    @cached_property
    def mfe_monomer(self):
        self._calculate_thermodynamics()
        return self.mfe_monomer

    @cached_property
    def mfe_homodimer(self):
        self._calculate_thermodynamics()
        return self.mfe_homodimer

    def _calculate_thermodynamics(self):
        self.mfe_monomer, self.mfe_homodimer = viennarna.dna_secondary_structure(
            self.seq
        )

    def __str__(self):
        return str(sequence.get_seq(self.seq))

    def __repr__(self):
        tm = self.__dict__.get("tm", "None")
        mfe_monomer = self.__dict__.get("mfe_monomer", "None")
        mfe_homodimer = self.__dict__.get("mfe_homodimer", "None")
        primer3 = "{...}" if self.primer3 else "None"
        return f"Primer(overhang={self.overhang}, binding={self.binding}, tm={tm}, mfe_monomer={mfe_monomer}, mfe_homodimer={mfe_monomer}, primer3={primer3})"


# class PrimerPair(NamedTuple):
#     primer1: Primer
#     primer2: Primer
#     ta: Optional[float]
#     mfe_heterodimer: Optional[float]

#     @property
#     def min_length(self):
#         return min(len(self.primer1.seq), len(self.primer2.seq))

#     @property
#     def max_length(self):
#         return max(len(self.primer1.seq), len(self.primer2.seq))


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
    overhang="",
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
            primer_seq = overhang + binding_seq
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


def get_primer_overhang(template, primer):
    sites = find_primer_binding_site(template, primer, return_sequences=True)
    if len(sites) != 1:
        raise ValueError(
            f"expecting a unique primer binding site, instead got {len(sites)}"
        )
    sense = sites[0][0]  # TODO: NamedTuple-ize
    score = sites[0][2]
    primer_overlap = sites[0][4]
    binding = primer_overlap[len(primer_overlap) - score :]
    overhang = primer[: len(primer) - len(primer_binding)]
    return overhang, binding, sense


def replace_primer_overhang(template, primer, new_overhang=""):
    overhang, binding, sense = get_primer_overhang(template, primer)
    if sense == -1:
        new_overhang = reverse_complement(new_overhang)
    new_primer = new_overhang + binding
    return new_primer


def format_primer_orientation(template, primer, format="f"):
    uppercase_initial = format[0].isupper()
    uppercase_all = format[-1].isupper()
    format = format.lower()
    sites = find_primer_binding_site(template, primer)
    if len(sites) != 1:
        raise ValueError(
            f"expecting a unique primer binding site, instead got {len(sites)}"
        )
    sense = sites[0][0]
    if format in ("r", "f"):
        if sense == -1:
            s = "r"
        else:
            s = "f"
    elif format in ("rev", "fwd"):
        if sense == -1:
            s = "rev"
        else:
            s = "fwd"
    elif format in ("reverse", "forward"):
        if sense == -1:
            s = "reverse"
        else:
            s = "forward"
    else:
        raise ValueError(f"unknown format string {format}")
    if uppercase_all:
        return s.upper()
    else:
        if uppercase_initial:
            return s.capitalize()
        else:
            return s
