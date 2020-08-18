import io
from Bio import SeqIO


def read_sequence(data, filename=None):
    buf = io.StringIO(data)
    if filename is None or filename.endswith("gb") or filename.endswith("gbk"):
        seq = SeqIO.read(buf, "genbank")
    elif filename.endswith("dna"):
        seq = SeqIO.read(buf, "snapgene")
    else:
        raise ValueError(f"unknown extension: {filename}")
    return seq
