import io
import requests
from Bio import SeqIO


def get_genbank(url):
    res = requests.get(url)
    buf = io.StringIO(res.content.decode("utf8"))
    seq = SeqIO.read(buf, "genbank")
    return seq
