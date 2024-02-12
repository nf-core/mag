#!/usr/bin/env python

## Originally written by Sabrina Krakau and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import sys
import argparse
import os.path
import pandas as pd
import csv
import gzip
import statistics

from Bio import SeqIO


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b",
        "--bins",
        required=True,
        nargs="+",
        metavar="FILE",
        help="Bins: FASTA containing all contigs.",
    )
    parser.add_argument(
        "-d",
        "--depths",
        required=True,
        metavar="FILE",
        help="(Compressed) TSV file containing contig depths for each sample: contigName, contigLen, totalAvgDepth, sample1_avgDepth, sample1_var [, sample2_avgDepth, sample2_var, ...].",
    )
    parser.add_argument(
        "-a", "--assembler", required=True, type=str, help="Assembler name."
    )
    parser.add_argument(
        "-i", "--id", required=True, type=str, help="Sample or group id."
    )
    parser.add_argument(
        "-m", "--binner", required=True, type=str, help="Binning method."
    )
    return parser.parse_args(args)


# Processing contig depths for each binner again, i.e. not the most efficient way, but ok


def main(args=None):
    args = parse_args(args)

    # load contig depths for all samples into dict (could use pandas as well)
    sample_names = []
    dict_contig_depths = {}
    with gzip.open(args.depths, "rt") as infile:
        reader = csv.reader(infile, delimiter="\t")
        # process header
        header = next(reader)
        for sample in range(int((len(header) - 3) / 2)):
            col_name = header[3 + 2 * sample]
            # retrieve sample name: "<assembler>-<id>-<other sample_name>.bam"
            sample_name = col_name[len(args.assembler) + 1 + len(args.id) + 1 : -4]
            sample_names.append(sample_name)
        # process contig depths
        for row in reader:
            contig_depths = []
            for sample in range(int((len(row) - 3) / 2)):
                contig_depths.append(float(row[3 + 2 * sample]))
            dict_contig_depths[str(row[0])] = contig_depths

    # Initialize output files
    n_samples = len(sample_names)
    with open(
        args.assembler + "-" + args.binner + "-" + args.id + "-binDepths.tsv", "w"
    ) as outfile:
        print("bin", "\t".join(sample_names), sep="\t", file=outfile)

    # for each bin, access contig depths and compute mean bin depth (for all samples)
    for file in args.bins:
        all_depths = [[] for i in range(n_samples)]

        if file.endswith(".gz"):
            with gzip.open(file, "rt") as infile:
                for rec in SeqIO.parse(infile, "fasta"):
                    contig_depths = dict_contig_depths[rec.id]
                    for sample in range(n_samples):
                        all_depths[sample].append(contig_depths[sample])
        else:
            with open(file, "rt") as infile:
                for rec in SeqIO.parse(infile, "fasta"):
                    contig_depths = dict_contig_depths[rec.id]
                    for sample in range(n_samples):
                        all_depths[sample].append(contig_depths[sample])

        binname = os.path.basename(file)
        with open(
            args.assembler + "-" + args.binner + "-" + args.id + "-binDepths.tsv", "a"
        ) as outfile:
            print(
                binname,
                "\t".join(
                    str(statistics.median(sample_depths))
                    for sample_depths in all_depths
                ),
                sep="\t",
                file=outfile,
            )


if __name__ == "__main__":
    sys.exit(main())
