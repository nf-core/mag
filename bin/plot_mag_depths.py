#!/usr/bin/env python

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
    parser.add_argument('-d', '--bin_depths'   , required=True, metavar='FILE'          , help="Bin depths file in TSV format (for one assembly): bin, sample1_depth, sample2_depth, ....")
    parser.add_argument('-g', '--groups'       , required=True, metavar='FILE'          , help="File in TSV format containing group information for samples: sample, group")
    parser.add_argument('-o', "--out"          , required=True, metavar='FILE', type=str, help="Output file.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # load data
    df     = pd.read_csv(args.bin_depths, sep='\t', index_col=0)
    groups = pd.read_csv(args.groups, sep='\t', index_col=0, names=['sample', 'group'])

    # add pseudo-abundances (sample-wise? dependent on lib-size)
    pseudo_cov = 0.1 * df[df > 0].min().min()
    df.replace(0, pseudo_cov, inplace=True)
    print(df)
    # compute centered log-ratios
    # divide df by sample-wise geometric means
    gmeans = stats.gmean(df, axis=0)                # apply on axis=0: 'index'
    df = np.log(df.div(gmeans, axis='columns'))     # divide column-wise (axis=1|'columns'), take natural logorithm
    # NOTE NaNs should not occur
    print(df)
    df.index.name='MAGs'
    df.columns.name='Samples'

    # prepare colors for group information
    color_map= dict(zip(groups['group'].unique(), sns.color_palette(n_colors=len(groups['group'].unique()))))

    # plot
    plt.figure()
    # TODO row_cluster=False ? add colors for tax. classification at certain level, e.g. phylum (if adding clustering based on phylogenetic information)
    # TODO yticklabels=False ? if number of bins > 20 ?
    sns.clustermap(df, row_cluster=True, cmap="vlag", col_colors=groups.group.map(color_map))
    plt.savefig(args.out)


if __name__ == "__main__":
    sys.exit(main())

