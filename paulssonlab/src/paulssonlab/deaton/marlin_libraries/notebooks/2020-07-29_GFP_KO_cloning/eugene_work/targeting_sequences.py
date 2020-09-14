# Add 10 bp before and after the sequence of interest
# Adjust the orientation of the sequence of interest to be from left to right
# Must be a fasta file
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import sys
import csv
import random
import copy
from Bio import pairwise2
import distance
import numpy as np
from Genome import PAM, bad_seq, comp_PAM, bad_target

# Input is a list of strings
def bad_seed(mut_list):
    updated_list = []
    size = len(mut_list)
    # Moving all bad_seeds into list called bad_seed_library
    # Change!
    with open("/Users/pang.eugene/CRISPR_og_lib/bad_seed_list.csv", "r") as bad_seed:
        reader = csv.DictReader(bad_seed)
        bad_seed_library = []
        for row in reader:
            bad_seed_library.append(row["seeds"])
    for x in range(size):
        seq = mut_list[x]
        if seq[15:20].upper not in bad_seed_library:
            updated_list.append(seq)
    if updated_list == []:
        return None
    return updated_list


# Counts number of mutations on a seq
def num_mut(seq):
    nucleotides = ["a", "t", "g", "c"]
    count = 0
    for x in seq:
        if x in nucleotides:
            count += 1
    return count


# Identifying targeting sequences with overlap and get rid of all of them
def overlap(library):
    updated_list = copy.deepcopy(library)
    size = len(updated_list)
    for x in range(size - 1):
        for y in range(x + 1, 20):
            target1 = library[x]
            target2 = library[y]
            mismatch = distance.hamming(target1, target2)
            alignment = 20 - mismatch
            if alignment > 8 and library[x] in updated_list:
                updated_list.remove(library[x])
                if library[y] in updated_list:
                    updated_list.remove(library[y])
                break
    return updated_list


# Turns a sequence of DNA into numbers
def str2int(string):
    code = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
    conv_str = np.array(list(map(lambda x: code[x], string)))
    return conv_str


# When I can figure out how to input l
# if len(sys.argv) == 2:
# record = SeqIO.read(sys.argv[1], "fasta")
# print(record)
# else:
# print("Error! Format: python target_sequences.py file.fasta")
# quit()
# Change!
record = SeqIO.read("/Users/pang.eugene/CRISPR_og_lib/gfpmut2.fasta", "fasta")
og_seq = record.seq
og_revc_seq = og_seq.reverse_complement()
seq = str(og_seq.upper())
revc_seq = str(og_revc_seq.upper())
library = []  # Reference library, will not be altered
bad_seed_library = []
count = 0
pam = ["AGG", "TGG", "GGG", "CGG"]
nucleotides = ["a", "t", "g", "c"]

# Identifying all of the targeting sequences 5' -> 3'
for x in range(len(seq)):
    if revc_seq[x : (x + 3)] in pam and ((x - 20) > 0):
        targeting = revc_seq[(x - 20) : (x)]
        library.append(targeting)
library = overlap(library)
size = len(library)

# Nested dictionary for each targeting with number of mismatches as the key
num_dict_list = []
for s in range(size):
    value = [None for i in range(11)]
    key = [
        "original",
        "one",
        "two",
        "three",
        "four",
        "five",
        "six",
        "seven",
        "eight",
        "nine",
        "ten",
    ]
    zipObj = zip(key, value)
    num_dict_list.append(dict(zipObj))
dictionary = dict(zip(library, num_dict_list))

# Fill in original only if it's not a bad seed
for target in library:
    result = bad_seed([target])
    dictionary[target]["original"] = result

# For 1 bp mismatch
for targeting in library:
    mut_list = []
    for x in range(20):
        for y in nucleotides:
            if y.upper != targeting[x]:
                mut_target = targeting[0:x] + y + targeting[x + 1 : 20]
                mut_list.append(mut_target)
    mut_list = bad_seed(mut_list)
    dictionary[targeting]["one"] = mut_list

# For 2+ bp mismatch will make 500 random and filter out bad seeds then select random 50
noobers = list(range(20))
for target in library:
    # number = number of mismatches
    for number in range(2, 11):
        mut_list = []
        for x in range(500):
            rand_noobers = random.sample(noobers, len(noobers))
            seq = target
            for y in range(number):
                n = rand_noobers[y]
                p = random.choice(nucleotides)
                while seq[n] == p.upper():
                    p = random.choice(nucleotides)
                seq = seq[0:n] + p + seq[n + 1 : 20]
            mut_list.append(seq)
        mut_list = bad_seed(mut_list)
        dictionary[target][key[number]] = mut_list

# Combine all lists for each target in library, preserves the index locations in indicies
# indicies order is library[0]key[0], library[0]key[1], etc
total = []
indicies = []
end = 0
for target in library:
    for x in key:
        start = end
        end = len(dictionary[target][x])
        indicies.append([start, end])
        total = total + dictionary[target][x]

# Finding bad targets in Genome.py
# Change!
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
size_bad_targets = len(bad_targets)

# Checking all of the potential targeting strands to the bad_targets
# Remove all potential targets with 10 matches or more
# 1st change all sequences into numbers for easier comparison
int_bad_targets = np.array(list(map(str2int, bad_targets)), dtype="uint8")
int_targets = np.array(list(map(str2int, total)), dtype="uint8")
# Reshaping for broadcast operation (N1, L) -> (1, N1, L)
int_bad_targets_b = np.array(int_bad_targets[np.newaxis, :, :])
# Reshaping for broadcast operation (N2, L) -> (N2, 1, L)
int_targets_b = np.array(int_targets[:, np.newaxis, :])
# Broadcast comparison (N2, N1, L)
bool_arr = int_targets_b == int_bad_targets_b
# Summing over L (N2, N1), calculate number of mismatches
match_arr = np.sum(bool_arr, axis=2)
good_targets_indicies = []
for y in range(len(total)):
    # True (good) when the number of matches is less than 10
    row_bool = map(lambda x: x > 9, match_arr[y, :])
    # If none are False
    if np.any(row_bool) == False:
        good_targets_indicies.append(y)

# Extract from large list back into dictionary
