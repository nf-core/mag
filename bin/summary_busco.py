#!/usr/bin/env python

# USAGE: ./summary.busco.py -s <summaries> -f <failed_bins>

import re
import sys
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', "--summaries", nargs="+", metavar='FILE', help="List of Busco summary files.")
    parser.add_argument('-f', "--failed_bins", nargs="+", metavar='FILE', help="List of files containing bin name for which Busco analysis failed.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if not args.summaries and not args.failed_bins:
        sys.exit("Either --summaries or --failed_bins must be specified!")

    # "# Summarized benchmarking in BUSCO notation for file MEGAHIT-testset1.contigs.fa"
    # "	C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:148"

    regexes = [r"# Summarized benchmarking in BUSCO notation for file (\S+)", r"	C:(\S+)%\[S:",
            r"%\[S:(\S+)%,D:", r"%,D:(\S+)%\],F:", r"%\],F:(\S+)%,M:", r"%,M:(\S+)%,n:", r"%,n:(\S+)"]
    columns = ["GenomeBin", "%Complete", "%Complete and single-copy",
            "%Complete and duplicated", "%Fragmented", "%Missing", "Total number"]

    # Search each summary file using its regex
    print("\t".join(columns))
    if args.summaries:
        for FILE in args.summaries:
            with open(FILE) as x:
                results = []
                TEXT = x.read()
                for REGEX in regexes:
                    match = re.search(REGEX, TEXT)
                    if match:
                        results.append(match.group(1))
                print("\t".join(results))

    # Add entries for bins with failed analysis
    if args.failed_bins:
        for FILE in args.failed_bins:
            with open(FILE) as x:
                failed_bin = x.readline().split('\n')[0]
                print(failed_bin, "0.0%", "0.0%", "0.0%", "0.0%", "100.0%", "NA", sep='\t')


if __name__ == "__main__":
    sys.exit(main())