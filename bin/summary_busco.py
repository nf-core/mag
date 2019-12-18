#!/usr/bin/env python

# USAGE: ./summary.busco.py *.txt

import re
from sys import argv

# "# Summarized benchmarking in BUSCO notation for file MEGAHIT-testset1.contigs.fa"
# "	C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:148"

regexes = [r"# Summarized benchmarking in BUSCO notation for file (\S+)", r"	C:(\S+)%\[S:",
           r"%\[S:(\S+)%,D:", r"%,D:(\S+)%\],F:", r"%\],F:(\S+)%,M:", r"%,M:(\S+)%,n:", r"%,n:(\S+)"]
columns = ["GenomeBin", "%Complete", "%Complete and single-copy",
           "%Complete and duplicated", "%Fragmented", "%Missing", "Total number"]

# Search each file using its regex
print("\t".join(columns))
for FILE in argv[1:]:
    with open(FILE) as x:
        results = []
        TEXT = x.read()
        for REGEX in regexes:
            match = re.search(REGEX, TEXT)
            if match:
                results.append(match.group(1))
        print("\t".join(results))
