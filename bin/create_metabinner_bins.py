#!/usr/bin/env python

## Originally written by Hesham Almessady (@HeshamAlmessady) and Adrian Fritz (@AlphaSquad)
# in https://github.com/hzi-bifo/mag and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import gzip
import io
import os
import sys

from Bio import SeqIO


def main():
    # Argument parsing
    if len(sys.argv) != 6:
        print(
            "Usage: python create_metabinner_bins.py <binning_file> <fasta_file> <output_path> <prefix> <length_threshold>"
        )
        sys.exit(1)

    binning = sys.argv[1]
    fasta = sys.argv[2]
    path = sys.argv[3]
    prefix = sys.argv[4]
    length = int(sys.argv[5])

    root = os.path.dirname(os.path.normpath(path)) or "."
    os.makedirs(path, exist_ok=True)

    metabinner_bins = {}
    with open(binning, "r") as b:
        for line in b:
            contig, bin = line.strip().split("\t")
            metabinner_bins[contig] = bin

    handles = {}

    def get_handle(dest_dir, fname):
        key = os.path.join(dest_dir, fname)
        if key not in handles:
            raw = open(key + ".gz", "wb")
            gz = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
            handles[key] = (io.TextIOWrapper(gz, encoding="utf-8"), raw)
        return handles[key][0]

    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record) < length:
                out = get_handle(root, prefix + ".tooShort.fa")
            elif record.id not in metabinner_bins:
                out = get_handle(root, prefix + ".unbinned.fa")
            else:
                out = get_handle(path, prefix + "." + metabinner_bins[record.id] + ".fa")
            SeqIO.write(record, out, "fasta")

    for text, raw in handles.values():
        text.close()
        raw.close()


if __name__ == "__main__":
    main()
