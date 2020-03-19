import io
import requests
from Bio import SeqIO


def get_genbank(url):
    res = requests.get(url)
    gb = res.content
    buf = io.StringIO(gb.decode("utf8"))
    dna = list(SeqIO.parse(buf, "genbank"))
    return dna
