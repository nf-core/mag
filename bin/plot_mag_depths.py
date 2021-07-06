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
    parser.add_argument('-o', "--out"          , required=True, metavar='FILE', type=str, help="Output file.")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # load data
    df=pd.read_csv(args.bin_depths, sep='\t', index_col=0)
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

    # plot
    plt.figure()
    # TODO add colors for group information for samples
    # TODO col_cluster=False ?
    # TODO yticklabels=False ? if number of bins > 20 ?
    # TODO add colors for tax. classification at certain level, e.g. phylum (if adding clustering based on phylogenetic information)
    sns.clustermap(df, row_cluster=False, cmap="vlag")
    plt.savefig(args.out)


if __name__ == "__main__":
    sys.exit(main())

