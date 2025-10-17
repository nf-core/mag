#!/usr/bin/env python

# Originally written by Sabrina Krakau and released under the MIT license.
# See git repository (https://github.com/nf-core/mag) for full license text.

import sys
import argparse
import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--bin_depths",
        required=True,
        metavar="FILE",
        help="Bin depths file in TSV format (for one assembly and binning method): bin, sample1_depth, sample2_depth, ....",
    )
    parser.add_argument(
        "-g",
        "--groups",
        required=True,
        metavar="FILE",
        help="File in TSV format containing group information for samples: sample, group",
    )
    parser.add_argument(
        "-o", "--out", required=True, metavar="FILE", type=str, help="Output file."
    )
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # load data
    df = pd.read_csv(args.bin_depths, sep="\t", index_col=0)
    groups = pd.read_csv(args.groups, sep="\t", index_col=0, names=["sample", "group"])

    # add pseudo-abundances (sample-wise? dependent on lib-size)
    pseudo_cov = 0.1 * df[df > 0].min().min()
    df.replace(0, pseudo_cov, inplace=True)
    # compute centered log-ratios
    # divide df by sample-wise geometric means
    gmeans = stats.gmean(df, axis=0)  # apply on axis=0: 'index'
    df = np.log(
        df.div(gmeans, axis="columns")
    )  # divide column-wise (axis=1|'columns'), take natural logorithm
    df.index.name = "MAGs"
    df.columns.name = "Samples"

    # prepare colors for group information
    color_map = dict(
        zip(
            groups["group"].unique(),
            sns.color_palette(n_colors=len(groups["group"].unique())),
        )
    )

    # plot
    plt.figure()
    bin_labels = True
    if len(df) > 30:
        bin_labels = False
    sns.clustermap(
        df,
        row_cluster=True,
        yticklabels=bin_labels,
        cmap="vlag",
        center=0,
        col_colors=groups.group.map(color_map),
        figsize=(6, 6),
    )
    plt.savefig(args.out)


if __name__ == "__main__":
    sys.exit(main())
