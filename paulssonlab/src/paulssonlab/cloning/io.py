import io
import os

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from paulssonlab.cloning.sequence import DsSeqRecord

TYPE_TO_MIMETYPE = {
    SeqRecord: "chemical/seq-na-genbank",
    Seq: "chemical/seq-na-fasta",
    str: "text/plain",
}

MIMETYPE_TO_EXTENSION = {
    "chemical/seq-na-genbank": ("gbk", "gb"),
    "chemical/seq-na-fasta": ("fasta",),
    "application/vnd.dna": ("dna",),
    # SEE: https://wiki.debian.org/DebianMedMIME
    "application/vnd.appliedbiosystems.abif": ("ab1",),
    "text/plain": ("txt",),
    None: (None,),  # bypass extension handling for None values
}
EXTENSION_TO_MIMETYPE = {}
for mimetype, exts in MIMETYPE_TO_EXTENSION.items():
    for ext in exts:
        EXTENSION_TO_MIMETYPE[ext] = mimetype


def value_to_bytes(value):
    if isinstance(value, (Seq, SeqRecord)):
        return value.format("genbank").encode()
    elif isinstance(value, str):
        return value.encode()
    else:
        raise ValueError(f"cannot convert {value.__class__} to bytes")


def _postprocess_sequence(seq):
    molecule_type = seq.annotations.get("molecule_type")
    if molecule_type is None:
        return seq
    elif molecule_type in ("ds-DNA", "DNA"):
        return DsSeqRecord(seq)
    else:
        raise NotImplementedError(
            f"cannot import genbank molecule_type: '{molecule_type}'"
        )


def bytes_to_value(bytes_, mimetype):
    if mimetype == "chemical/seq-na-genbank":
        buf = io.StringIO(bytes_.decode())
        return _postprocess_sequence(SeqIO.read(buf, "genbank"))
    elif mimetype == "application/vnd.dna":
        buf = io.BytesIO(bytes_)
        return _postprocess_sequence(SeqIO.read(buf, "snapgene"))
    elif mimetype == "chemical/seq-na-fasta":
        buf = io.StringIO(bytes_.decode())
        return SeqIO.read(buf, "fasta")
    elif mimetype == "application/vnd.appliedbiosystems.abif":
        buf = io.BytesIO(bytes_)
        return SeqIO.read(buf, "abi")
    elif mimetype == "text/plain":
        return bytes_.decode()
    raise ValueError(f"cannot convert mimetype {mimetype} from bytes")


def value_to_mimetype(value):
    if value is None:
        return None
    for type_, mimetype in TYPE_TO_MIMETYPE.items():
        if isinstance(value, type_):
            return mimetype
    raise ValueError(f"cannot find mimetype for {value.__class__}")


def filename_to_mimetype(filename):
    return EXTENSION_TO_MIMETYPE[os.path.splitext(filename)[1][1:].lower()]


def value_to_extension(value):
    return MIMETYPE_TO_EXTENSION[value_to_mimetype(value)][0]


def read_file(filename):
    with open(filename, "rb") as f:
        bytes_ = f.read()
    mimetype = filename_to_mimetype(filename)
    return bytes_to_value(bytes_, mimetype)


def read_http(url):
    mimetype = filename_to_mimetype(url)
    bytes_ = requests.get(url).content
    return bytes_to_value(bytes_, mimetype)
