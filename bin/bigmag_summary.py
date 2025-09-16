#!/usr/bin/env python

import pandas as pd
import re
import argparse
import sys
import warnings

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--summary", metavar="FILE", help="Pipeline summary file.")
    parser.add_argument("-g", "--gunc_summary", metavar="FILE", help="GUNC summary file.")
    parser.add_argument("-a", "--alt_summary", metavar="FILE", help="BUSCO or CheckM2 summary file.")
    parser.add_argument(
        "-t", "--binqc_tool", help="Bin QC tool used", choices=["busco", "checkm", "checkm2"]
    )

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
        and not args.alt_summary
    ):
        sys.exit(
            "No summary specified! "
            "Please specify the pipeline summary, the GUNC summary and BUSCO or CheckM2 summary."
        )

    df_summary = pd.read_csv(args.summary, sep='\t', index_col=0)
    for i in range(len(df_summary["bin"])):
        name = df_summary["bin"][i]
        name = re.sub(r'\.(fa|fasta)(\..*)?$', '', name)
        df_summary.at[i,"bin"] = name
        df_summary = df_summary.sort_values(by='bin')
        df_summary["bin"] = df_summary["bin"].astype(str)
    
    df_gunc = pd.read_csv(args.gunc_summary, sep='\t')
    df_gunc["genome"] = df_gunc["genome"].astype(str)
    df_gunc = df_gunc.sort_values(by='genome')

    df_alt = pd.read_csv(args.alt_summary, sep='\t')

    column_names = ['genome']

    if args.binqc_tool == "busco":
        df_alt["Name"] = df_alt["Name"].astype(str)
        df_alt = df_alt.sort_values(by='Name')
        column_names.append("Name")

    elif args.binqc_tool == "checkm" or args.binqc_tool == "checkm2":
        for i in range(len(df_alt["Input_file"])):
            name = df_alt["Input_file"][i]
            name = re.sub(r'\.(fa|fasta)(\..*)?$', '', name)
            df_alt.at[i,"Input_file"] = name
            df_alt = df_alt.sort_values(by='Input_file')
            df_alt["Input_file"] = df_alt["Input_file"].astype(str)
        column_names.append("Input_file")

    df_list = [df_gunc, df_alt]
    for i in range(len(df_list)):
        df_summary = pd.merge(df_summary, df_list[i], left_on='bin', right_on=column_names[i], how='left')
    
    df_summary.rename(columns={'bin': 'Bin'}, inplace=True)
    
    columns_to_remove = ['Name', "genome", 'Input_file', 'Assembly', 'Bin Id']
    for column in df_summary.columns:
        if column in columns_to_remove:
            df_summary = df_summary.drop(columns=column)

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

