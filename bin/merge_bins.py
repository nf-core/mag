#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse

from shutil import copyfile, copyfileobj

from Bio import SeqIO
from Bio.SeqUtils import GC


class Bin(object):
    """a genome bin"""
    def __init__(self, id, compl, contam):
        super(Bin, self).__init__()
        self.id = id
        self.compl = compl
        self.contam = contam
        self.read_cov = float()
        self.phylo = []


class PossibleMerge(object):
    """A possible merger between two bins. Corresponds to one line of
    merger.tsv"""
    def __init__(self, line):
        self.bin_1 = Bin(line[0], line[2], line[3])
        self.bin_2 = Bin(line[1], line[4], line[5])
        self.delta_compl = float(line[6])
        self.delta_contam = float(line[7])
        self.delta_merger = float(line[8])  # ???
        self.merged_compl = float(line[9])
        self.merged_contam = float(line[10])

    def add_profile_info(self, profile_txt, mean):
        for profile in self.parse_profile(profile_txt):
            cov = (float(profile[2]) * mean) / (float(profile[1]) * 1e6)
            cov = round(cov, 2)
            if profile[0] == self.bin_1.id and self.bin_1.read_cov == 0.0:
                self.bin_1.read_cov = cov
            elif profile[0] == self.bin_2.id and self.bin_2.read_cov == 0.0:
                self.bin_2.read_cov = cov

    def add_tree_info(self, tree_qa):
        for phylogeny in self.parse_tree(tree_qa):
            if phylogeny[0] == self.bin_1.id and self.bin_1.phylo == []:
                self.bin_1.phylo = phylogeny[3:]
            elif phylogeny[0] == self.bin_2.id and self.bin_2.phylo == []:
                self.bin_2.phylo = phylogeny[3:]

    @staticmethod
    def parse_profile(profile_txt):
        with open(profile_txt, "r") as p:
            p.readline()
            header = p.readline()
            p.readline()
            for l in p:
                sl = l.split()
                if len(sl) < 2:
                    continue
                yield(sl[:3])  # bin_id, bin size, n mapped reads

    @staticmethod
    def parse_tree(tree_qa):
        with open(tree_qa, "r") as t:
            t.readline()
            header = t.readline()
            t.readline()
            for l in t:
                sl = l.split()
                if len(sl) < 2:
                    continue
                yield(sum([i.split(";") for i in sl], []))


def get_mean(length_file):
    with open(length_file, "r") as f:
        total = [int(line.split()[0]) for line in f]
    mean = sum(total) / len(total)
    return(mean)


def parse_merge_file(merger_tsv):
    with open(merger_tsv, "r") as m:
        m.readline()  # discard the header
        for l in m:
            line = l.split()
            pm = PossibleMerge(line)
            yield(pm)


def concatenate(file_list, output):
    """Concatenate files together
    Args:
        file_list (list): the list of input files (can be a generator)
        output (string): the output file name
    """
    try:
        out_file = open(output, "wb")
    except (IOError, OSError) as e:
        print("Failed to open output file: %s" % e)
        sys.exit(1)

    with out_file:
        for file_name in file_list:
            if file_name is not None:
                with open(file_name, "rb") as f:
                    copyfile(f, out_file)


def merge(args):
    all_bins = os.listdir(args.bins)
    n_merged = 0

    # merge bins according to criteria in
    # https://doi.org/10.1038/s41564-017-0012-7
    mean_read_length = get_mean(args.length)
    for pm in parse_merge_file(args.merger):
        pm.add_profile_info(args.profile, mean_read_length)
        pm.add_tree_info(args.tree)
        abs_delta_cov = \
            max(pm.bin_1.read_cov, pm.bin_2.read_cov) /\
            min(pm.bin_1.read_cov, pm.bin_2.read_cov)
        clades = []
        for p1, p2 in zip(pm.bin_1.phylo, pm.bin_2.phylo):
            clades.append((p1, p2))

        if pm.delta_contam <= 5 and \
           pm.merged_contam <= 20 and \
           pm.delta_compl >= 10 and \
           abs_delta_cov <= 0.25 and \
           clades[-1][0] == clades[-1][1]:

            # then read the two bins
            bin_1_p = "%s/%s.fa" % (args.bins, pm.bin_1.id)
            bin_2_p = "%s/%s.fa" % (args.bins, pm.bin_2.id)
            with open(bin_1_p, "r") as bin_1, open(bin_2_, "r") as bin_2:
                concat_seq_1 = concat_seq_2 = ""
                for record in SeqIO.parse(bin_1, "fasta"):
                    concat_seq_1 += record.seq
                gc_bin_1 = GC(concat_seq_1)
                for record in SeqIO.parse(bin_2, "fasta"):
                    concat_seq_2 += record.seq
                gc_bin_2 = GC(concat_seq_2)

            # then merge if similar GC content
            if abs(gc_bin_1 - gc_bin_2) <= 3:
                output_name = "%s/merged_%i.fa" % (args.outdir, n_merged)
                concatenate([bin_1_p, bin_2_p], output_name)
                all_bins.remove(bin_1_p)
                all_bins.remove(bin_2_p)
                n_merged += 1

        else:
            continue
    # now puts the single bins in the same output directory
    for bin in all_bins:
        if bin.endswith(".fa") or bin.endswith(".fasta"):
            path = "%s/%s" % (args.bins, bin)
            output_name = "%s/%s" % (args.outdir, bin)
            copyfile(path, output_name)


def main():
    parser = argparse.ArgumentParser(
        prog="merge_bins.py",
        usage="parse merger.tsv from checkm and merge bins"
    )
    parser.add_argument(
        "--profile",
        metavar="profile.txt",
        help="output of checkm profile"
    )
    parser.add_argument(
        "--tree",
        metavar="tree_qa.txt",
        help="output of checkm tree_qa"
    )
    parser.add_argument(
        "--length",
        metavar="read_length.txt",
        help="a file containing read length distribution"
    )
    parser.add_argument(
        "--merger",
        metavar="merger.tsv",
        help="merger.tsv produced by checkm merge"
    )
    parser.add_argument(
        "bins",
        metavar="bins/",
        help="the bin directory (bins must end with .fa or .fasta)"
    )
    parser.add_argument(
        "outdir",
        metavar="outdir/",
        help="output directory (must exist)"
    )
    parser.set_defaults(func=merge)
    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError as e:
        parser.print_help()
        raise


if __name__ == "__main__":
    main()
