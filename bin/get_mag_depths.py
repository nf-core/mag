#!/usr/bin/env python

## Originally written by Sabrina Krakau and updated by Diego Alvarez. Released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import argparse
import csv
import gzip
import os.path
import statistics
import sys

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
        help=(
            "Produces (compressed) TSV file containing contig depths for each sample: contigName, contigLen, "
            + "totalAvgDepth, sample1_avgDepth, sample1_var [, sample2_avgDepth, sample2_var, ...]."
        ),
    )
    parser.add_argument("-a", "--assembler", required=True, type=str, help="Assembler name.")
    parser.add_argument("-i", "--id", required=True, type=str, help="Sample or group id.")
    parser.add_argument("-m", "--binner", required=True, type=str, help="Binning method.")
    return parser.parse_args(args)


# Processing contig depths for each binner again, i.e. not the most efficient way, but ok


def main(args=None):
    args = parse_args(args)

    # load contig depths for all samples into dict (could use pandas as well)
    sample_names = []
    dict_contig_depths = {}
    with gzip.open(args.depths, "rt") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        # process header to extract sample names from column names
        depth_columns = []
        for col_name in reader.fieldnames[3::2]:  # Every 2nd column starting from index 3
            # retrieve sample name: "<assembler>-<id>-<other sample_name>.bam"
            sample_name = col_name[len(args.assembler) + 1 + len(args.id) + 1 : -4]
            sample_names.append(sample_name)
            depth_columns.append(col_name)
        # process contig depths
        for row in reader:
            contig_depths = {}
            for sample_name, col_name in zip(sample_names, depth_columns):
                contig_depths[sample_name] = float(row[col_name])
            dict_contig_depths[row[reader.fieldnames[0]]] = contig_depths

    sample_names = sorted(sample_names)

    with open(f"{args.assembler}-{args.binner}-{args.id}-binDepths.tsv", "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["bin", *sample_names])

        # for each bin, access contig depths and compute mean bin depth (for all samples)
        for file in args.bins:
            all_depths = {sample: [] for sample in sample_names}

            if file.endswith(".gz"):
                with gzip.open(file, "rt") as infile:
                    for rec in SeqIO.parse(infile, "fasta"):
                        contig_depths = dict_contig_depths[rec.id]
                        for sample in sample_names:
                            all_depths[sample].append(contig_depths[sample])
            else:
                with open(file, "rt") as infile:
                    for rec in SeqIO.parse(infile, "fasta"):
                        contig_depths = dict_contig_depths[rec.id]
                        for sample in sample_names:
                            all_depths[sample].append(contig_depths[sample])

            binname = os.path.basename(file)

            writer.writerow(
                [binname, *[statistics.median(all_depths[sample]) for sample in sample_names]],
            )


if __name__ == "__main__":
    sys.exit(main())
