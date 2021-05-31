#!/usr/bin/env python

import sys
import argparse
import os.path
import pandas as pd
import csv
import gzip

from Bio import SeqIO


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bins'       , required=True, nargs="+", metavar='FILE'                             , help="Bins: FASTA containing all contigs.")
    parser.add_argument('-d', '--depths'     , required=True           , metavar='FILE'                             , help="(Compressed) TSV file containing contig depths for each sample (within group): contigName, contigLen, totalAvgDepth, sample1_avgDepth, sample1_var [, sample2_avgDepth, sample2_var, ...].")
    parser.add_argument('-o', "--out"        , required=True           , metavar='FILE', type=argparse.FileType('w'), help="Output file containing depth for each bin of assembly.")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    # load all contig depths into dict, roughly estimated the memory consuption should be feasible even for extremely large assemblies (?)
    # contig : [depth1, depth2, ...]

    sample_names = []
    dict_contig_depths = {}
    with gzip.open(args.depths, "rt") as infile:
        reader = csv.reader(infile, delimiter = "\t")
        # process header
        header = next(reader)
        for sample in range(int((len(header)-3)/2)):
            sample_names.append(header[3+2*sample])
        # process rest
        for row in reader:
            contig_depths = []
            for sample in range(int((len(row)-3)/2)):
                contig_depths.append(float(row[3+2*sample]))
            dict_contig_depths[str(row[0])] = contig_depths

    n_samples = len(sample_names)
    print("bin", '\t'.join(sample_names), file=args.out)
    for file in args.bins:
        mean_depths = [0] * n_samples
        with open(file, "rt") as infile:
            c = 0
            for rec in SeqIO.parse(infile,'fasta'):
                print("rec: ", rec.id)
                depths = dict_contig_depths[rec.id]
                for sample in range(n_samples):
                    mean_depths[sample] += depths[sample]
                c += 1
        for sample in range(n_samples):
            mean_depths[sample] = mean_depths[sample]/float(c)
        print(os.path.basename(file), '\t'.join(str(d) for d in mean_depths), file=args.out)


if __name__ == "__main__":
    sys.exit(main())
