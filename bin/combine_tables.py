#!/usr/bin/env python

import sys
import argparse
import os.path
import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--depths_summary", required=True, metavar="FILE", help="Bin depths summary file.")
    parser.add_argument("-b", "--busco_summary", metavar="FILE", help="BUSCO summary file.")
    parser.add_argument("-c", "--checkm_summary", metavar="FILE", help="CheckM summary file.")
    parser.add_argument("-q", "--quast_summary", metavar="FILE", help="QUAST BINS summary file.")
    parser.add_argument("-g", "--gtdbtk_summary", metavar="FILE", help="GTDB-Tk summary file.")

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing final summary.",
    )
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if not args.busco_summary and not args.checkm_summary and not args.quast_summary and not args.gtdbtk_summary:
        sys.exit("No summary specified! Please specify at least BUSCO, CheckM or QUAST summary.")

    # GTDB-Tk can only be run in combination with BUSCO or CheckM
    if args.gtdbtk_summary and not (args.busco_summary or args.checkm_summary):
        sys.exit("Invalid parameter combination: GTDB-TK summary specified, but no BUSCO or CheckM summary!")

    # handle bin depths
    results = pd.read_csv(args.depths_summary, sep="\t")
    results.columns = ["Depth " + str(col) if col != "bin" else col for col in results.columns]
    bins = results["bin"].sort_values().reset_index(drop=True)

    if args.busco_summary:
        busco_results = pd.read_csv(args.busco_summary, sep="\t")
        if not bins.equals(busco_results["GenomeBin"].sort_values().reset_index(drop=True)):
            sys.exit("Bins in BUSCO summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, busco_results, left_on="bin", right_on="GenomeBin", how="outer"
        )  # assuming depths for all bins are given

    if args.checkm_summary:
        checkm_results = pd.read_csv(args.checkm_summary, sep="\t")
        checkm_results["Bin Id"] = checkm_results["Bin Id"] + ".fa"
        if not bins.equals(checkm_results["Bin Id"].sort_values().reset_index(drop=True)):
            sys.exit("Bins in CheckM summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, checkm_results, left_on="bin", right_on="Bin Id", how="outer"
        )  # assuming depths for all bins are given
        results["Bin Id"] = results["Bin Id"].str.removesuffix(".fa")

    if args.quast_summary:
        quast_results = pd.read_csv(args.quast_summary, sep="\t")
        if not bins.equals(quast_results["Assembly"].sort_values().reset_index(drop=True)):
            sys.exit("Bins in QUAST summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, quast_results, left_on="bin", right_on="Assembly", how="outer"
        )  # assuming depths for all bins are given

    if args.gtdbtk_summary:
        gtdbtk_results = pd.read_csv(args.gtdbtk_summary, sep="\t")
        if not bins.equals(gtdbtk_results["user_genome"].sort_values().reset_index(drop=True)):
            sys.exit("Bins in GTDB-Tk summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, gtdbtk_results, left_on="bin", right_on="user_genome", how="outer"
        )  # assuming depths for all bins are given

    results.to_csv(args.out, sep="\t")


if __name__ == "__main__":
    sys.exit(main())
