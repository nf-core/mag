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

# Read file
with open(input_file) as f:
    fasta_sequences = SeqIO.parse(f,'fasta')

# make table
    df = pd.DataFrame(columns=['id','seq','length'])
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        length = len(sequence)
        df = df.append({"id":name, "seq":sequence, "length":length}, ignore_index = True)

#sort table by sequence length
df.sort_values(by=['length'], ascending=False)

# Save sequences
remaining = []
pooled = []
for index, row in df.iterrows():
    # each sequence above threshold into one file
    if row['length'] >= length_threshold and index+1 <= max_sequences:
        print("write "+out_base+"."+str(index+1)+".fa")
        out = (SeqRecord(Seq(row['seq'], generic_dna), id = row['id']))
        SeqIO.write(out, out_base+"."+str(index+1)+".fa", "fasta")
    # contigs to retain
    elif row['length'] >= min_length_to_retain_contig:
        pooled.append(SeqRecord(Seq(row['seq'], generic_dna), id = row['id']))
    # remaining sequences
    else:
        remaining.append(SeqRecord(Seq(row['seq'], generic_dna), id = row['id']))

print("write "+out_base+".pooled.fa")
SeqIO.write(pooled, out_base+".pooled.fa", "fasta")
print("write "+out_base+".remaining.fa")
SeqIO.write(remaining, out_base+".remaining.fa", "fasta")
