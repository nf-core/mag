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
    parser.add_argument('-b', '--bins'         , required=True, nargs="+", metavar='FILE'                             , help="Bins: FASTA containing all contigs.")
    parser.add_argument('-d', '--depths'       , required=True           , metavar='FILE'                             , help="(Compressed) TSV file containing contig depths for each sample: contigName, contigLen, totalAvgDepth, sample1_avgDepth, sample1_var [, sample2_avgDepth, sample2_var, ...].")
    parser.add_argument('-a', '--assembly_name', required=True                           , type=str                   , help="Assembly name.")
    parser.add_argument('-o', "--out"          , required=True           , metavar='FILE', type=argparse.FileType('w'), help="Output file containing depth for each bin.")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    # load contig depths for all samples into dict (could use pandas as well)
    sample_names = []
    dict_contig_depths = {}
    with gzip.open(args.depths, "rt") as infile:
        reader = csv.reader(infile, delimiter = "\t")
        # process header
        header = next(reader)
        for sample in range(int((len(header)-3)/2)):
            col_name = header[3+2*sample]
            # retrieve sample name: "<assembly_name>-<sample_name>.bam"
            sample_name = col_name[len(args.assembly_name)+1:-4]
            sample_names.append(sample_name)
        # process contig depths
        for row in reader:
            contig_depths = []
            for sample in range(int((len(row)-3)/2)):
                contig_depths.append(float(row[3+2*sample]))
            dict_contig_depths[str(row[0])] = contig_depths

    n_samples = len(sample_names)
    # for each bin, access contig depths and compute mean bin depth (for all samples)
    print("bin", '\t'.join(sample_names), sep='\t', file=args.out)
    for file in args.bins:
        mean_depths = [0] * n_samples
        with open(file, "rt") as infile:
            c = 0
            for rec in SeqIO.parse(infile,'fasta'):
                depths = dict_contig_depths[rec.id]
                for sample in range(n_samples):
                    mean_depths[sample] += depths[sample]
                c += 1
        for sample in range(n_samples):
            mean_depths[sample] = mean_depths[sample]/float(c)
        print(os.path.basename(file), '\t'.join(str(d) for d in mean_depths), sep='\t', file=args.out)


if __name__ == "__main__":
    sys.exit(main())
