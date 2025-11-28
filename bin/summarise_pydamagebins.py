#!/usr/bin/env python3

## summarise_pydamagebins.py
## Originally written by James Fellows Yates and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import os
import sys
import argparse
import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--contig_to_bin_map",
        required=True,
        metavar="FILE",
        help="Input file of tab-separated list of <bin_id>\t<contig_id> mappings. No header.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        metavar="FILE",
        help="Path to output TSV file with all statistic values summarised with as a median value per bin.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        action=argparse.BooleanOptionalAction,
        help="Print to console more logging information",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument(
        "pydamage_reports", nargs="+", help="List of pydamage report files for contigs."
    )

    return parser.parse_args(args)


def main(args=None):
    parser = argparse.ArgumentParser(description="summarise_pydamage.py")
    args = parse_args(args)

    if args.output is not None:
        if os.path.abspath(args.output) == os.path.abspath(args.input):
            sys.exit(
                "[summarise_pydamage.py] ERROR: Input and output files names must be different."
            )

    if args.verbose:
        print(
            "[summarise_pydamagebins.py] PROCESSING: Loading " + args.contig_to_bin_map
        )

    contig_to_bin_map = pd.read_csv(
        args.contig_to_bin_map, sep="\t", names=["bin_id", "contig_id"]
    )

    pydamage_reports_dfs = []
    for f in args.pydamage_reports:
        if args.verbose:
            print("[summarise_pydamagebins.py] PROCESSING: Loading " + f)
        assembly_id = os.path.basename(f).replace("_pydamage_results.csv", "")
        pydamage_df = pd.read_csv(f, sep=",")
        pydamage_df["assembly_id"] = assembly_id
        pydamage_df.insert(0, "assembly_id", pydamage_df.pop("assembly_id"))
        pydamage_reports_dfs.append(pydamage_df)

    pydamage_reports_all = pd.concat(pydamage_reports_dfs)

    ## tODO:
    ##    - Split assembly_id intwo two columns
    ##    - Same for contig_to_bin_map
    ##    - Merge on contig_id / reference
    ##    - Summarise

    ## Testing only
    pydamage_reports_all.to_csv("pydamage_bin_summary.tsv", sep="\t", index=False)
