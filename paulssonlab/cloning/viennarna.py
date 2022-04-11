import os
from threading import Lock
from paulssonlab.util.threading import synchronized
import RNA

PARAMETERS = {
    "dna_mathews2004": f"{os.environ['CONDA_PREFIX']}/share/ViennaRNA/dna_mathews2004.par"
}

viennarna_lock = Lock()
current_parameter = None


def configure_parameters(param_file):
    global current_parameter
    if current_parameter == param_file:
        return
    else:
        RNA.read_parameter_file(param_file)
        current_parameter = param_file


# FROM: https://github.com/ViennaRNA/ViennaRNA/issues/64
@synchronized(viennarna_lock)
def dna_secondary_structure(seq):
    configure_parameters(PARAMETERS["dna_mathews2004"])
    md = RNA.md()
    fc = RNA.fold_compound(seq, md)
    (_, mfe_monomer) = fc.mfe()
    fc_dimer = RNA.fold_compound(f"{seq}&{seq}", md)
    (_, mfe_homodimer) = fc_dimer.mfe()
    return mfe_monomer, mfe_homodimer


@synchronized(viennarna_lock)
def dna_heterodimer(seq1, seq2):
    configure_parameters(PARAMETERS["dna_mathews2004"])
    md = RNA.md()
    fc_dimer = RNA.fold_compound(f"{seq1}&{seq2}", md)
    (_, mfe_heterodimer) = fc_dimer.mfe()
    return mfe_heterodimer
