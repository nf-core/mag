#!/usr/bin/env python

## Originally written by Jeferyd Yepes and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import argparse
import re
import sys

import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--summary", metavar="FILE", help="Pipeline summary file.")
    parser.add_argument("-g", "--gunc_summary", metavar="FILE", help="GUNC summary file.")

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="FILE",
        type=argparse.FileType("w"),
        help="Output file containing final bigmag summary.",
    )
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    if (
        not args.summary
        and not args.gunc_summary
    ):
        sys.exit(
            "No summary specified! "
            "Please specify the pipeline summary and the GUNC summary."
        )

    df_summary = pd.read_csv(args.summary, sep='\t')
    df_summary.columns = df_summary.columns.str.replace(r'(_busco|_checkm2|_checkm|_gtdbtk|_gunc|_quast)$', '', regex=True)
    for i in range(len(df_summary["bin"])):
        name = df_summary["bin"][i]
        name = re.sub(r'\.(fa|fasta)(\..*)?$', '', name)
        df_summary.at[i,"bin"] = name
        df_summary = df_summary.sort_values(by='bin')
        df_summary["bin"] = df_summary["bin"].astype(str)

    df_gunc = pd.read_csv(args.gunc_summary, sep='\t')
    df_gunc["genome"] = df_gunc["genome"].astype(str)
    df_gunc = df_gunc.sort_values(by='genome')

    df_summary = pd.merge(df_summary, df_gunc, left_on='bin', right_on='genome', how='left')

    df_summary.rename(columns={'bin': 'Bin'}, inplace=True)
    columns_to_remove = ['Name', "genome", 'Input_file', 'Assembly', 'Bin Id']
    df_summary = df_summary.drop(columns=columns_to_remove, errors="ignore")

    df_summary['sample'] = None
    for f in range(len(df_summary["Bin"])):
        match = re.search(r'^.*?-.*?-(.*)$', df_summary["Bin"][f])
        if match:
            name = match.group(1)
            name = re.sub(r'\.(unbinned|noclass)(\..*)?$', '', name)
            name = re.sub(r'\.\d+(\.[^.]+)?$', '', name)
            df_summary.at[f,"sample"] = name

    df_summary.to_csv(args.out, sep="\t", index=True)

if __name__ == "__main__":
    sys.exit(main())
