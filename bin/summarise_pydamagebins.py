#!/usr/bin/env python3

## summarise_pydamagebins.py
## Originally written by James Fellows Yates and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import os
import sys
import argparse
import pandas as pd
from pathlib import Path


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
        "pydamage_reports_input",
        nargs="+",
        help="List of pydamage report files for contigs.",
    )

    return parser.parse_args(args)


def main(args=None):
    argparse.ArgumentParser(description="summarise_pydamage.py")
    args = parse_args(args)

    if args.verbose:
        print(
            "[summarise_pydamagebins.py] PROCESSING: Loading " + args.contig_to_bin_map
        )

    ## Load contig to bin map from Nextflow
    contig_to_bin_map = pd.read_csv(
        args.contig_to_bin_map,
        sep="\t",
    )
    contig_to_bin_map.rename(columns={"contig_id": "reference"}, inplace=True)
    contig_to_bin_map["assembly_id"] = contig_to_bin_map["assembly_id"].astype(str)

    ## Repair contig names to match
    ## Some tools remove everything after first space, so pyDamage has truncated headers vs. raw contigs
    ## We strip this information from the raw contigs to meet at the simplest/lowest common denominator
    contig_to_bin_map["reference"] = contig_to_bin_map.apply(
        lambda row: (row["reference"].split(" ")[0]),
        axis=1,
    )

    ## Clean up pydamage reports
    pydamage_reports_dfs = []

    for report in args.pydamage_reports_input:
        if args.verbose:
            print("[summarise_pydamagebins.py] PROCESSING: Loading " + report)
        assembly_id = os.path.basename(report).replace("_pydamage_results.csv", "")
        pydamage_df = pd.read_csv(report, sep=",")
        pydamage_df["assembly_id"] = assembly_id
        pydamage_df.insert(0, "assembly_id", pydamage_df.pop("assembly_id"))
        pydamage_reports_dfs.append(pydamage_df)

    pydamage_reports_all = pd.concat(pydamage_reports_dfs)
    pydamage_reports_all["assembly_id"] = pydamage_reports_all["assembly_id"].astype(
        str
    )

    ## Merge contig_to_bin_map with pydamage reports
    if args.verbose:
        print(
            "[summarise_pydamagebins.py] MERGING: contig to bin map with pydamage reports"
        )

    pydamage_contig_to_bin = pd.merge(
        left=contig_to_bin_map,
        right=pydamage_reports_all,
        on=["assembly_id", "reference"],
    ).sort_values(by=["bin_id", "assembly_id", "binner", "reference"])

    ## Group by bin_id, save per-bin collection, and then summarise over median
    if args.verbose:
        print("[summarise_pydamagebins.py] GROUPING: by bin and calculating median")

    pydamage_bin_summary = pydamage_contig_to_bin.groupby(["bin_id"], as_index=False)

    for name, group in pydamage_bin_summary:
        filename = Path(name).stem + "_pydamage_bin_results.tsv"
        group.to_csv(filename, sep="\t", index=False)

    pydamage_bin_summary_median = pydamage_bin_summary.median(
        numeric_only=True
    ).sort_values(by=["bin_id"])

    ## Testing only
    pydamage_bin_summary_median.to_csv(
        "pydamage_bins_summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    sys.exit(main())
