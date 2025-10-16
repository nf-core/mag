#!/usr/bin/env python

## Originally written by Daniel Straub and Sabrina Krakau and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import argparse
import sys
import warnings

import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--depths_summary",
        required=True,
        metavar="FILE",
        help="Bin depths summary file.",
    )
    parser.add_argument("-b", "--binqc_summary", metavar="FILE", help="BUSCO summary file.")
    parser.add_argument("-q", "--quast_summary", metavar="FILE", help="QUAST BINS summary file.")
    parser.add_argument("-g", "--gtdbtk_summary", metavar="FILE", help="GTDB-Tk summary file.")
    parser.add_argument("-a", "--cat_summary", metavar="FILE", help="CAT table file.")
    parser.add_argument(
        "-t", "--binqc_tool", help="Bin QC tool used", choices=["busco", "checkm", "checkm2"]
    )

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing final summary.",
    )
    return parser.parse_args(args)


def parse_cat_table(cat_table):
    """Parse CAT table.

    CAT table is trickier to parse than the other tables, because it has a variable number of columns,
    depending on the number of ranks that are reported for the taxonomic assignation of each contig.
    Therefore, we first parse the header to get the column names, and then parse the table, to get the
    maximum number of columns. Then, we merge the columns containing the ranks into a single column.

    Args:
        cat_table (str): Path to CAT table

    Returns:
        pd.DataFrame: parse CAT table
    """
    with open(cat_table, "r") as f:
        next(f)  # skip header
        maxcol = 0
        for line in f:
            maxcol = max(maxcol, len(line.split("\t")))

    header = [
        "bin",
        "classification",
        "reason",
        "lineage",
        "lineage scores",
        "full lineage names",
    ]

    df = pd.read_table(
        cat_table,
        names=header + [f"rank_{i}" for i in range(maxcol - len(header))],
        on_bad_lines="warn",
        header=None,
        skiprows=1,
    )
    # merge all rank columns into a single column
    df["CAT_rank"] = (
        df.filter(regex="rank_\d+").apply(lambda x: ";".join(x.dropna()), axis=1).str.lstrip()
    )
    # remove rank_* columns
    df.drop(df.filter(regex="rank_\d+").columns, axis=1, inplace=True)

    return df


def main(args=None):
    args = parse_args(args)

    if (
        not args.binqc_summary
        and not args.quast_summary
        and not args.gtdbtk_summary
    ):
        sys.exit(
            "No summary specified! "
            "Please specify at least BUSCO, CheckM, CheckM2 or QUAST summary."
        )

    # GTDB-Tk can only be run in combination with BUSCO, CheckM or CheckM2
    if args.gtdbtk_summary and not args.binqc_summary:
        sys.exit(
            "Invalid parameter combination: "
            "GTDB-TK summary specified, but no BUSCO, CheckM or CheckM2 summary!"
        )

    # handle bin depths
    results = pd.read_csv(args.depths_summary, sep="\t")
    results.columns = ["Depth " + str(col) if col != "bin" else col for col in results.columns]
    bins = results["bin"].sort_values().reset_index(drop=True)

    if args.binqc_summary and args.binqc_tool == "busco":
        busco_results = pd.read_csv(args.binqc_summary, sep="\t")
        busco_bins = set(busco_results["Input_file"])

        if set(bins) != busco_bins and len(busco_bins.intersection(set(bins))) > 0:
            warnings.warn("Bins in BUSCO summary do not match bins in bin depths summary")
        elif len(busco_bins.intersection(set(bins))) == 0:
            sys.exit("Bins in BUSCO summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, busco_results, left_on="bin", right_on="Input_file", how="outer"
        )  # assuming depths for all bins are given

    if args.binqc_summary and args.binqc_tool == "checkm":
        use_columns = [
            "Bin Id",
            "Marker lineage",
            "# genomes",
            "# markers",
            "# marker sets",
            "Completeness",
            "Contamination",
            "Strain heterogeneity",
            "Coding density",
            "Translation table",
            "# predicted genes",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5+",
        ]
        checkm_results = pd.read_csv(args.binqc_summary, usecols=use_columns, sep="\t")
        checkm_results["Bin Id"] = checkm_results["Bin Id"] + ".fa"
        if not set(checkm_results["Bin Id"]).issubset(set(bins)):
            sys.exit("Bins in CheckM summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, checkm_results, left_on="bin", right_on="Bin Id", how="outer"
        )  # assuming depths for all bins are given
        results["Bin Id"] = results["Bin Id"].str.removesuffix(".fa")

    if args.binqc_summary and args.binqc_tool == "checkm2":
        use_columns = [
            "Name",
            "Completeness",
            "Contamination",
            "Completeness_Model_Used",
            "Coding_Density",
            "Translation_Table_Used",
            "Total_Coding_Sequences",
        ]
        checkm2_results = pd.read_csv(args.binqc_summary, usecols=use_columns, sep="\t")
        checkm2_results["Name"] = checkm2_results["Name"] + ".fa"
        if not set(checkm2_results["Name"]).issubset(set(bins)):
            sys.exit("Bins in CheckM2 summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, checkm2_results, left_on="bin", right_on="Name", how="outer"
        )  # assuming depths for all bins are given
        results["Name"] = results["Name"].str.removesuffix(".fa")

    if args.quast_summary:
        quast_results = pd.read_csv(args.quast_summary, sep="\t")
        if not bins.equals(quast_results["Assembly"].sort_values().reset_index(drop=True)):
            sys.exit("Bins in QUAST summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, quast_results, left_on="bin", right_on="Assembly", how="outer"
        )  # assuming depths for all bins are given

    if args.gtdbtk_summary:
        gtdbtk_results = pd.read_csv(args.gtdbtk_summary, sep="\t")
        if len(set(gtdbtk_results["user_genome"].to_list()).difference(set(bins))) > 0:
            sys.exit("Bins in GTDB-Tk summary do not match bins in bin depths summary!")
        results = pd.merge(
            results, gtdbtk_results, left_on="bin", right_on="user_genome", how="outer"
        )  # assuming depths for all bins are given

    if args.cat_summary:
        cat_results = parse_cat_table(args.cat_summary)
        if len(set(cat_results["bin"].to_list()).difference(set(bins))) > 0:
            sys.exit("Bins in CAT summary do not match bins in bin depths summary!")
        results = pd.merge(
            results,
            cat_results[["bin", "CAT_rank"]],
            left_on="bin",
            right_on="bin",
            how="outer",
        )

    results.to_csv(args.out, sep="\t")


if __name__ == "__main__":
    sys.exit(main())
