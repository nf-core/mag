#!/usr/bin/env python

import sys
import argparse
import os.path
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', "--busco_summary",                    metavar='FILE',                               help="BUSCO summary file.")
    parser.add_argument('-q', "--quast_summary",                    metavar='FILE',                               help="QUAST BINS summary file.")
    parser.add_argument('-g', "--gtdbtk_summary",                   metavar='FILE',                               help="GTDB-Tk summary file.")

    parser.add_argument('-o',  "--out",               required=True, metavar='FILE', type=argparse.FileType('w'), help="Output file containing final summary.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if not args.busco_summary and not args.quast_summary and not args.gtdbtk_summary:
        sys.exit("No summary specified! Please specify at least two summaries.")

    # At least two summaries must be specified for merging
    # Since GTDB-Tk can only be run in combination with BUSCO, test if (BUSCO && GTDB-TK) || (BUSCO && QUAST)
    if args.gtdbtk_summary and not args.busco_summary:
        sys.exit("Invalid parameter combination: GTDB-TK summary specified, but no BUSCO summary!")
    if args.quast_summary and not args.busco_summary:
        sys.exit("At least two summaries must be specified, but only QUAST summary was given! Please specify BUSCO summary.")

    results = pd.DataFrame()
    bins    = pd.Series()

    if args.busco_summary:
        results = pd.read_csv(args.busco_summary, sep="\t")
        bins    = results['GenomeBin'].sort_values().reset_index(drop=True)

    if args.quast_summary:
        quast_results = pd.read_csv(args.quast_summary, sep="\t")
        if bins.empty:
            bins = quast_results['Assembly'].sort_values().reset_index(drop=True)
        elif not bins.equals(quast_results['Assembly'].sort_values().reset_index(drop=True)):
            sys.exit("Bins in QUAST summary do not match bins in BUSCO summary!")
        results = pd.merge(results, quast_results, left_on="GenomeBin", right_on="Assembly", how='outer')

    if args.gtdbtk_summary:
        gtdbtk_results = pd.read_csv(args.gtdbtk_summary, sep="\t")
        if not bins.equals(gtdbtk_results['user_genome'].sort_values().reset_index(drop=True)):
            sys.exit("Bins in GTDB-Tk summary do not match bins in BUSCO summary!")   # GTDB-Tk can currently anyway only run in combination with BUSCO
        results = pd.merge(results, gtdbtk_results, left_on="GenomeBin", right_on="user_genome", how='outer')   # again assuming BUSCO summary must be given

    # Print to stdout
    results.to_csv(args.out, sep='\t')


if __name__ == "__main__":
    sys.exit(main())
