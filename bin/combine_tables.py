#!/usr/bin/env python

#USAGE: ./combine_tables.py <BUSCO_table> <QUAST_table>

import pandas as pd
from sys import stdout
from sys import argv

# Read files
file1 = pd.read_csv(argv[1], sep="\t")
file2 = pd.read_csv(argv[2], sep="\t")

# Merge files
result = pd.merge(file1, file2, left_on="GenomeBin", right_on="Assembly", how='outer')

# Print to stdout
result.to_csv(stdout, sep='\t')
