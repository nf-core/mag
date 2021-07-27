#!/usr/bin/env python

#USAGE: ./combine_tables.py <*.unbinned.fa> <length threshold> <maximal number of sequences> <length threshold to retain contigs>

import pandas as pd
from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import os

# Input
input_file = argv[1]
length_threshold = int(argv[2])
max_sequences = int(argv[3])
min_length_to_retain_contig = int(argv[4])

# Base name for file output
out_base = (os.path.splitext(input_file)[0])

# Data structures to separate and store sequences
df_above_threshold = pd.DataFrame(columns=['id','seq','length'])
pooled             = []
remaining          = []

# Read file
with open(input_file) as f:
    fasta_sequences = SeqIO.parse(f,'fasta')

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        length = len(sequence)

        # store each sequence above threshold together with its length into df
        if length >= length_threshold:
            df_above_threshold = df_above_threshold.append({"id":name, "seq":sequence, "length":length}, ignore_index = True)
        # contigs to retain and pool
        elif length >= min_length_to_retain_contig:
            pooled.append(SeqRecord(Seq(sequence, generic_dna), id = name))
        # remaining sequences
        else:
            remaining.append(SeqRecord(Seq(sequence, generic_dna), id = name))

# Sort sequences above threshold by length
df_above_threshold.sort_values(by=['length'], ascending=False, inplace=True)
df_above_threshold.reset_index(drop=True, inplace=True)

# Write `max_sequences` longest sequences (above threshold) into separate files, add remainder to pooled
for index, row in df_above_threshold.iterrows():
    if index+1 <= max_sequences:
        print("write "+out_base+"."+str(index+1)+".fa")
        out = (SeqRecord(Seq(row['seq'], generic_dna), id = row['id']))
        SeqIO.write(out, out_base+"."+str(index+1)+".fa", "fasta")
    else:
        pooled.append(SeqRecord(Seq(row['seq'], generic_dna), id = row['id']))

print("write "+out_base+".pooled.fa")
SeqIO.write(pooled, out_base+".pooled.fa", "fasta")
print("write "+out_base+".remaining.fa")
SeqIO.write(remaining, out_base+".remaining.fa", "fasta")
