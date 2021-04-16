#!/usr/bin/env python

# USAGE: ./summary.busco.py -sd <summaries_domain> -ss <summaries_specific> -f <failed_bins>

import re
import sys
import argparse
import os.path
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--auto', default=False, action='store_true', help="BUSCO run in auto lineage selection mode.")
    parser.add_argument('-sd', "--summaries_domain", nargs="+", metavar='FILE', help="List of BUSCO summary files for domains.")
    parser.add_argument('-ss', "--summaries_specific", nargs="+", metavar='FILE', help="List of BUSCO summary files for specific lineages.")
    parser.add_argument('-f', "--failed_bins", nargs="+", metavar='FILE', help="List of files containing bin name for which BUSCO analysis failed.")
    parser.add_argument('-o', "--out", required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing final BUSCO summary.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if not args.summaries_domain and not args.summaries_specific and not args.failed_bins:
        sys.exit("Either --summaries_domain, --summaries_specific or --failed_bins must be specified!")

    # "# Summarized benchmarking in BUSCO notation for file /path/to/MEGAHIT-testset1.contigs.fa"
    # "	C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:148"

    regexes = [r"# Summarized benchmarking in BUSCO notation for file (\S+)", r"# The lineage dataset is: (\S+) \(", r"	C:(\S+)%\[S:",
            r"%\[S:(\S+)%,D:", r"%,D:(\S+)%\],F:", r"%\],F:(\S+)%,M:", r"%,M:(\S+)%,n:", r"%,n:(\S+)"]
    columns_domain = ["GenomeBin", \
                     "Domain", \
                     "%Complete (domain)", \
                     "%Complete and single-copy (domain)", \
                     "%Complete and duplicated (domain)", \
                     "%Fragmented (domain)", \
                     "%Missing (domain)", \
                     "Total number (domain)"]
    columns_specific = ["GenomeBin", \
                        "Specific lineage dataset", \
                        "%Complete (specific)", \
                        "%Complete and single-copy (specific)", \
                        "%Complete and duplicated (specific)", \
                        "%Fragmented (specific)", \
                        "%Missing (specific)", \
                        "Total number (specific)"]

    if args.auto:
        columns = ["GenomeBin", \
                        "Domain", \
                        "%Complete (domain)", \
                        "%Complete and single-copy (domain)", \
                        "%Complete and duplicated (domain)", \
                        "%Fragmented (domain)", \
                        "%Missing (domain)", \
                        "Total number (domain)", \
                        "Specific lineage dataset", \
                        "%Complete (specific)", \
                        "%Complete and single-copy (specific)", \
                        "%Complete and duplicated (specific)", \
                        "%Fragmented (specific)", \
                        "%Missing (specific)", \
                        "Total number (specific)"]
    else:
        columns = ["GenomeBin", \
                        "Specific lineage dataset", \
                        "%Complete (specific)", \
                        "%Complete and single-copy (specific)", \
                        "%Complete and duplicated (specific)", \
                        "%Fragmented (specific)", \
                        "%Missing (specific)", \
                        "Total number (specific)"]

    # Search each summary file using its regex
    results_domain = []
    if args.summaries_domain:
        for file in args.summaries_domain:
            with open(file) as infile:
                results = []
                text = infile.read()
                for index, regex in enumerate(regexes):
                    match = re.search(regex, text)
                    if match:
                        if index == 0:
                            results.append(os.path.basename(match.group(1)))
                        else:
                            results.append(match.group(1))
                results_domain.append(results)
    df_domain = pd.DataFrame(results_domain, columns=columns_domain)

    results_specific = []
    if args.summaries_specific:
        for file in args.summaries_specific:
            with open(file) as infile:
                results = []
                text = infile.read()
                for index, regex in enumerate(regexes):
                    match = re.search(regex, text)
                    if match:
                        if index == 0:
                            results.append(os.path.basename(match.group(1)))
                        else:
                            results.append(match.group(1))
                results_specific.append(results)
    df_specific = pd.DataFrame(results_specific, columns=columns_specific)

    # Add entries for bins with failed analysis (for domain and specific lineage where applicable)
    failed = []
    if args.failed_bins:
        for file in args.failed_bins:
            with open(file) as infile:
                failed_bin = infile.readline().split('\n')[0]
                if args.auto:
                    results = [failed_bin, "NA", "0.0%", "0.0%", "0.0%", "0.0%", "100.0%", "NA", "NA", "0.0%", "0.0%", "0.0%", "0.0%", "100.0%", "NA"]
                else:
                    results = [failed_bin, "NA", "0.0%", "0.0%", "0.0%", "0.0%", "100.0%", "NA"]
                failed.append(results)
    df_failed = pd.DataFrame(failed, columns=columns)

    # merge results
    if args.auto:
        df_final = df_domain\
            .merge(df_specific, on="GenomeBin", how='outer')\
            .append(df_failed)
    else:
        df_final = df_specific\
            .append(df_failed)

    df_final.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())