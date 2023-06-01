import copy
import csv
import random
import re
import sys

import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature


# Find location of all PAM sites on sequence
def PAM(genome):
    locations = re.finditer("GG", genome)
    indices = [m.start(0) for m in locations]
    return indices


# Input the SeqI0.read result to find all locations the targeting strand has negative outcomes
# Template_strand is 1 or -1 depending on top or bottom strand respectively
def bad_seq(record, template_strand):
    bad_locations = []
    regulatory = ["promoter", "regulatory", "terminator", "protein_bind", "rep_origin"]
    coding_regulatory = ["tRNA", "rRNA", "gene", "CDS", "ncRNA"]
    count = 0
    for feat in record.features:
        L = feat.location
        count += 1
        if feat.type in regulatory and (count % 2 != 0):
            start = int(L.start)
            end = int(L.end)
            i = [start, end]
            bad_locations.append(i)
        elif feat.type in coding_regulatory and (count % 2 != 0):
            # Template strand must be opposite of coding strand
            if L.strand == -1 and template_strand == 1:
                start = int(L.start)
                end = int(L.end)
                i = [start, end]
                bad_locations.append(i)
            elif L.strand == 1 and template_strand == -1:
                start = int(L.start)
                end = int(L.end)
                i = [start, end]
                bad_locations.append(i)
    return bad_locations


# Compares to see if any of the PAMS are in the bad sequence regions
def comp_PAM(results_PAM, results_bad_seq):
    bad_pam = []
    bad_seq_array = np.array(results_bad_seq)
    for pam in results_PAM:
        above_start_arr = pam >= bad_seq_array[:, 0]
        below_end_arr = pam <= bad_seq_array[:, 1]
        test_arr = above_start_arr * below_end_arr
        bad_pam_test = np.any(test_arr)
        if bad_pam_test == True:
            bad_pam.append(pam)
    return bad_pam


# Generates all the bad pam targeting sequences
def bad_target(bad_pam_locations, genome):
    library = []
    for x in bad_pam_locations:
        if (x - 20) > 0:
            targeting = genome[(x - 20) : (x)]
            library.append(targeting)
    return library


# record = SeqIO.read("/Users/pang.eugene/CRISPR_og_lib/gfpmut2.fasta", "fasta")
record = SeqIO.read(
    "/Users/pang.eugene/CRISPR_og_lib/CRISPRi_reference_genome.gb", "genbank"
)
seq = record.seq

# Find location of all the PAM sequences
rev_pam = PAM((str(seq.reverse_complement())).upper())
fwd_pam = PAM((str(seq)).upper())

# Find all of the locations that are off limits
rev_locations = bad_seq(record.reverse_complement(), -1)
fwd_locations = bad_seq(record, 1)

# See if any of the PAM sequences are in the off limits locations
bad_pam_fwd = comp_PAM(fwd_pam, fwd_locations)
bad_pam_rev = comp_PAM(rev_pam, rev_locations)

# Building potential bad targeting strand library
bad_target_fwd = bad_target(bad_pam_fwd, (str(seq)).upper())
bad_target_rev = bad_target(bad_pam_rev, (str(seq.reverse_complement())).upper())
bad_targets = bad_target_fwd + bad_target_rev
