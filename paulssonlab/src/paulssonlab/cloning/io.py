import io
from Bio import SeqIO
from paulssonlab.cloning.sequence import DsSeqRecord


def read_sequence(data, filename=None):
    buf = io.StringIO(data)
    if filename is None or filename.endswith("gb") or filename.endswith("gbk"):
        seq = SeqIO.read(buf, "genbank")
    elif filename.endswith("dna"):
        seq = SeqIO.read(buf, "snapgene")
    else:
        raise ValueError(f"unknown extension: {filename}")
    molecule_type = seq.annotations["molecule_type"]
    if molecule_type == "ds-DNA":
        seq = DsSeqRecord(seq)
    else:
        raise NotImplementedError(
            f"cannot import genbank molecule_type: '{molecule_type}'"
        )
    return seq


def gdrive_file(content):
    # TODO: read/write ab1
    if isinstance(content, (DsSeqRecord, SeqRecord)):
        content = content.format("genbank")
        mimetype = "chemical/seq-na-genbank"
    elif isinstance(content, Seq):
        content = content.format("fasta")
        mimetype = "chemical/seq-na-fasta"
    elif isinstance(content, str):
        mimetype = "text/plain"
    else:
        raise ValueError(f"cannot upload type {content.__class__} to Google Drive")
    file = {"content": content, "mimetype": mimetype}
    return file
