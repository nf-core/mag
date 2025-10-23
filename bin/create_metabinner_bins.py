#!/usr/bin/env python

## Originally written by Hesham Almessady (@HeshamAlmessady) and Adrian Fritz (@AlphaSquad) in https://github.com/hzi-bifo/mag and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import sys
import os
from Bio import SeqIO

binning = sys.argv[1]
fasta = sys.argv[2]
path = sys.argv[3]
prefix = sys.argv[4]
length = sys.argv[5]

try:
    os.mkdir(path)
except OSerror as error:
    pass

Metabinner_bins = {}
with open(binning, 'r') as b:
    for line in b:
        contig, bin = line.strip().split('\t')
        Metabinner_bins[contig] = bin

for outf in set(Metabinner_bins.values()):
    name = prefix + "." + outf + ".fa"
    open(os.path.join(path, name), 'w')

with open(fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if len(record) < int(length):
            f = prefix + ".tooShort.fa"
        elif record.id not in Metabinner_bins:
            f = prefix + ".unbinned.fa"
        else:
            f = prefix + "." + Metabinner_bins[record.id] + ".fa"
        with open(os.path.join(path, f), 'a') as out:
            SeqIO.write(record, out, "fasta")
