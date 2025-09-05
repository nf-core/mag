#!/usr/bin/env python3

## summarise_pydamage.py
## Originally written by James Fellows Yates and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

import os
import sys
import argparse
import pandas as pd


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        help="Input CSV file output of a pyDamage analyze command",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        metavar="FILE",
        help="Path to output TSV file with all statistic values summarised with as a median value. If not supplied, default output is input with _summarised appended to the file name.",
    )
    parser.add_argument(
        "-n",
        "--name",
        required=True,
        metavar="STRING",
        help="Sample name for appending to summarised result as ID",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        action=argparse.BooleanOptionalAction,
        help="Print to console more logging information",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.0.1")

    return parser.parse_args(args)


def main(args=None):
    parser = argparse.ArgumentParser(description="summarise_pydamage.py")
    args = parse_args(args)

    if args.output is not None:
        if os.path.abspath(args.output) == os.path.abspath(args.input):
            sys.exit(
                "[summarise_pydamage.py] ERROR: Input and output files names must be different."
            )

    if args.verbose:
        print("[summarise_pydamage.py] PROCESSING: Loading " + args.input)

    pydamage_raw = pd.read_csv(args.input, sep=",")

    if args.verbose:
        print("[summarise_pydamage.py] PROCESSING: appending sample name")

    pydamage_raw["id"] = os.path.splitext(args.input)[0]

    if args.verbose:
        print(
            "[summarise_pydamage.py] PROCESSING: cleaning up table, and calculating median"
        )

    pydamage_summarised = pydamage_raw.drop("reference", axis=1).groupby("id").median()

    if args.verbose:
        print("[summarise_pydamage.py] FINALISING: saving file")

    if args.output is not None:
        if os.path.splitext(args.output)[1] != ".tsv":
            outfile = os.path.abspath(args.output + ".tsv")
        else:
            outfile = os.path.abspath(args.output)
    else:
        outfile = os.path.splitext(args.input)[0] + "_summarised.tsv"

    pydamage_summarised.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    sys.exit(main())
