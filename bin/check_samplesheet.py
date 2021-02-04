#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/design.csv


import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/mag samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO nf-core: Update the check_samplesheet function
def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample, group, short_reads_1, short_reads_2, long_reads
    # TODO ...
    """

    group_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 4
        # TODO nf-core: Update the column names for the input samplesheet
        HEADER = ["sample", "group", "short_reads_1", "short_reads_2", "long_reads"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample_name, group, sr1, sr2, lr = lspl[: len(HEADER)]
            if sample_name:
                if sample_name.find(" ") != -1:
                    print_error("Sample entry contains spaces!", "Line", line)
            else:
                print_error("Sample entry has not been specified!", "Line", line)
            ## Check group name entries
            if group:
                if group.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [sr1, sr2, lr]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, sr1, sr2, lr]
            if sr1 and sr2:
                sample_info = ["0", sr1, sr2, lr]
            elif sr1 and not sr2 and lr:
                print_error("Invalid combination of single-end short reads and long reads provided! SPAdes does not support single-end data and thus hybrid assembly cannot be performed.", "Line", line)
            elif sr1 and not sr2:
                sample_info = ["1", sr1, sr2, lr]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)
            ## Create sample mapping dictionary = {group: {sample_name : single_end, sr1, sr2 }}
            if group not in group_run_dict:
                group_run_dict[group] = {}
            if sample_name in group_run_dict.values():
                print_error("Samplesheet contains duplicate sample IDs!", "Line", line)
            else:
                group_run_dict[group][sample_name] = sample_info


    ## Write validated samplesheet with appropriate columns
    if len(group_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:

            fout.write(",".join(["sample", "group", "single_end", "short_reads_1", "short_reads_2", "long_reads"]) + "\n")
            for group in sorted(group_run_dict.keys()):

                for sample in sorted(group_run_dict[group].keys()):

                    ## Write to file
                    sample_info = group_run_dict[group][sample]
                    fout.write(",".join([sample] + [group] + sample_info) + "\n")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
