#!/usr/bin/env python

## Originally written by Hesham Almessady (@HeshamAlmessady) and Adrian Fritz (@AlphaSquad) in https://github.com/hzi-bifo/mag and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import sys
import os
from Bio import SeqIO

def main():
    # Argument parsing
    if len(sys.argv) != 6:
        print("Usage: python create_metabinner_bins.py <binning_file> <fasta_file> <output_path> <prefix> <length_threshold>")
        sys.exit(1)

    binning = sys.argv[1]
    fasta = sys.argv[2]
    path = sys.argv[3]
    prefix = sys.argv[4]
    length = int(sys.argv[5])

    # Create output directory if it doesn't exist
    os.makedirs(path, exist_ok=True)

    # Load binning data into a dictionary
    Metabinner_bins = {}
    with open(binning, 'r') as b:
        for line in b:
            contig, bin = line.strip().split('\t')
            Metabinner_bins[contig] = bin

    # Process the input fasta file
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record) < length:
                f = prefix + ".tooShort.fa"
            elif record.id not in Metabinner_bins:
                f = prefix + ".unbinned.fa"
            else:
                f = prefix + "." + Metabinner_bins[record.id] + ".fa"
            with open(os.path.join(path, f), 'a') as out:
                SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    main()
