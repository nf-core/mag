#!/usr/bin/env python

# Originally written by Sabrina Krakau and released under the MIT license.
# See git repository (https://github.com/nf-core/mag) for full license text.

import re
import sys
import argparse
import os.path
import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-x",
        "--extension",
        required=True,
        type=str,
        help="File extension passed to GTDB-TK and substracted by GTDB-Tk from bin names in results files.",
    )
    parser.add_argument(
        "-s",
        "--summaries",
        nargs="+",
        metavar="FILE",
        help="List of GTDB-tk summary files.",
    )
    parser.add_argument(
        "-fi",
        "--filtered_bins",
        nargs="+",
        metavar="FILE",
        help="List of files containing names of bins which where filtered out during GTDB-tk analysis.",
    )
    parser.add_argument(
        "-fa",
        "--failed_bins",
        nargs="+",
        metavar="FILE",
        help="List of files containing bin names for which GTDB-tk analysis failed.",
    )
    parser.add_argument(
        "-d",
        "--qc_discarded_bins",
        nargs="+",
        metavar="FILE",
        type=str,
        help="List of files containing names of bins which were discarded based on BUSCO metrics.",
    )

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing final GTDB-tk summary.",
    )
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if (
        not args.summaries
        and not args.filtered_bins
        and not args.failed_bins
        and not args.qc_discarded_bins
    ):
        sys.exit(
            "Either --summaries, --filtered_bins, --failed_bins or --qc_discarded_bins must be specified!"
        )

    columns = [
        "user_genome",
        "classification",
        "closest_genome_reference",
        "closest_genome_reference_radius",
        "closest_genome_taxonomy",
        "closest_genome_ani",
        "closest_genome_af",
        "closest_placement_reference",
        "closest_placement_radius",
        "closest_placement_taxonomy",
        "closest_placement_ani",
        "closest_placement_af",
        "pplacer_taxonomy",
        "classification_method",
        "note",
        "other_related_references(genome_id,species_name,radius,ANI,AF)",
        "msa_percent",
        "translation_table",
        "red_value",
        "warnings",
    ]
    # Note: currently all columns included

    # For bins already discarded based on BUSCO QC metrics
    discarded = []
    if args.qc_discarded_bins:
        for bin_name in args.qc_discarded_bins:
            bin_results = [
                bin_name,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
                pd.NA,
            ]
            discarded.append(bin_results)

    df_final = pd.DataFrame(discarded, columns=columns)
    df_final.set_index("user_genome", inplace=True)

    # For bins with succesfull GTDB-tk classification
    if args.summaries:
        for file in args.summaries:
            df_summary = pd.read_csv(file, sep="\t")[columns]
            # add by GTDB-Tk substracted file extension again to bin names (at least until changed consistently in rest of pipeline)
            df_summary["user_genome"] = (
                df_summary["user_genome"].astype(str) + "." + args.extension
            )
            df_summary.set_index("user_genome", inplace=True)
            df_final = df_final.append(df_summary, verify_integrity=True)

    # For bins that were filtered out by GTDB-tk (e.g. due to insufficient number of AAs in MSA)
    filtered = []
    if args.filtered_bins:
        for file in args.filtered_bins:
            df = pd.read_csv(file, sep="\t", names=["bin_name", "reason"])
            for index, row in df.iterrows():
                bin_name = row["bin_name"]
                bin_results = [
                    bin_name,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                ]
                filtered.append(bin_results)

    df_filtered = pd.DataFrame(filtered, columns=columns)
    df_filtered["user_genome"] = (
        df_filtered["user_genome"].astype(str) + "." + args.extension
    )
    df_filtered.set_index("user_genome", inplace=True)
    df_final = df_final.append(df_filtered, verify_integrity=True)

    # For bins for which GTDB-tk classification failed
    failed = []
    if args.failed_bins:
        for file in args.failed_bins:
            df = pd.read_csv(file, sep="\t", names=["bin_name", "reason"])
            for index, row in df.iterrows():
                bin_name = row["bin_name"]
                bin_results = [
                    bin_name,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                    pd.NA,
                ]
                failed.append(bin_results)

    df_failed = pd.DataFrame(failed, columns=columns)
    df_failed["user_genome"] = (
        df_failed["user_genome"].astype(str) + "." + args.extension
    )
    df_failed.set_index("user_genome", inplace=True)
    df_final = df_final.append(df_failed, verify_integrity=True)

    # write output
    df_final.reset_index().rename(columns={"index": "user_genome"}).to_csv(
        args.out, sep="\t", index=False
    )


if __name__ == "__main__":
    sys.exit(main())
