#!/usr/bin/env python

#USAGE: ./summary.busco.py *.txt

import re
from sys import argv

#"# Summarized benchmarking in BUSCO notation for file MEGAHIT-testset1.contigs.fa"
#"	C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:148"

regexes = {
    'nf-core/mag': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'fastqc': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'fastp': ['v_fastp.txt', r"fastp (\S+)"],
    'megahit': ['v_megahit.txt', r"MEGAHIT v(\S+)"],
    'metabat': ['v_metabat.txt', r"version (\S+)"],
    'NanoPlot': ['v_nanoplot.txt', r"NanoPlot (\S+)"],
    'Filtlong': ['v_filtlong.txt', r"Filtlong v(\S+)"],
    'porechop': ['v_porechop.txt', r"(\S+)"],
    'NanoLyse': ['v_nanolyse.txt', r"NanoLyse (\S+)"],
    'SPAdes': ['v_spades.txt', r"SPAdes v(\S+)"],
    'BUSCO': ['v_busco.txt', r"BUSCO (\S+)"],
}

regexes = [r"# Summarized benchmarking in BUSCO notation for file (\S+)", r"	C:(\S+)%\[S:", r"%\[S:(\S+)%,D:", r"%,D:(\S+)%\],F:", r"%\],F:(\S+)%,M:", r"%,M:(\S+)%,n:", r"%,n:(\S+)"]
columns = ["GenomeBin","%Complete","%Complete and single-copy","%Complete and duplicated","%Fragmented","%Missing","Total number"]

# Search each file using its regex
print("\t".join(columns))
for FILE in argv[1:]:
    with open(FILE) as x:
        results = []
        TEXT = x.read()
        for REGEX in regexes:
            match = re.search(REGEX, TEXT)
            if match:
                results.append( match.group(1) )
        print("\t".join(results))
