import numpy as np
from io import StringIO
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import tempfile


def read_ab1(ab1_file, array_wrapper=np.array):
    seq = SeqIO.read(ab1_file, "abi")
    quality = seq.letter_annotations["phred_quality"]
    peak_calls = seq.annotations["abif_raw"]["PLOC2"]
    order = seq.annotations["abif_raw"]["FWO_1"]
    sequence = str(seq.seq)
    # Extract color data
    peaks = {}
    for i, (base, r) in enumerate(zip(order, [9, 10, 11, 12])):
        peaks[base] = array_wrapper(seq.annotations["abif_raw"]["DATA{}".format(r)])
    return {
        "name": seq.id,
        "sequence": sequence,
        "peaks": peaks,
        "peak_calls": array_wrapper(peak_calls),
        "quality": array_wrapper(quality),
    }


def write_fasta(seqs, handle):
    # records = [SeqRecord(Seq(seq, IUPAC.DNA), name) for name, seq in seqs.items()]
    records = [SeqRecord(Seq(seq), name) for name, seq in seqs.items()]
    return SeqIO.write(records, handle, "fasta")


def call_subprocess(args, input_data=None):
    if input_data is not None:
        stdin = subprocess.PIPE
    else:
        stdin = None
    proc = subprocess.Popen(
        args, stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if input_data is not None:
        input_bytes = input_data.encode("utf-8")
    else:
        input_bytes = None
    stdout, stderr = proc.communicate(input_bytes)
    if stderr:
        raise Exception('got stderr: "{}"'.format(stderr))
    return stdout


def align_clustalomega(seqs):
    input_data = StringIO()
    write_fasta(seqs, input_data)
    output_data = call_subprocess(
        ["clustalo", "--infile=-", "--seqtype=DNA", "--outfmt=clustal"],
        input_data.getvalue(),
    )
    # return list(SeqIO.parse(StringIO(output_data.decode('utf-8')), 'clustal'))#, IUPAC.DNA)
    return next(
        AlignIO.parse(StringIO(output_data.decode("utf-8")), "clustal")
    )  # , IUPAC.DNA)


def align_mafft(seqs):
    with tempfile.NamedTemporaryFile() as input_handle:
        write_fasta(seqs, input_handle.name)  # TODO: write directly using file handle
        output_data = call_subprocess(
            ["mafft", "--quiet", "--adjustdirection", input_handle.name]
        )
        return next(
            AlignIO.parse(StringIO(output_data.decode("utf-8")), "fasta")
        )  # , IUPAC.DNA)


def trim_to_ref(msa, ref_name="ref"):
    ary = np.array([list(str(s.seq)) for s in msa if s.description != ref_name])
    mask = (ary == "-").sum(axis=0) >= len(msa) - 1
    idx_start = (~mask).argmax()
    idx_stop = (~mask)[::-1].argmax()
    return msa[:, idx_start:-idx_stop]
