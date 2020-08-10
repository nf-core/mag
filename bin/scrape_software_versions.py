#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/mag': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'fastp': ['v_fastp.txt', r"fastp (\S+)"],
    'MEGAHIT': ['v_megahit.txt', r"MEGAHIT v(\S+)"],
    'Metabat': ['v_metabat.txt', r"version (\S+)"],
    'NanoPlot': ['v_nanoplot.txt', r"NanoPlot (\S+)"],
    'Filtlong': ['v_filtlong.txt', r"Filtlong v(\S+)"],
    'Porechop': ['v_porechop.txt', r"(\S+)"],
    'NanoLyse': ['v_nanolyse.txt', r"NanoLyse (\S+)"],
    'SPAdes': ['v_spades.txt', r"SPAdes v(\S+)"],
    'BUSCO': ['v_busco.txt', r"BUSCO (\S+)"],
    'Centrifuge': ['v_centrifuge.txt', r"centrifuge-class version (\S+)"],
    'Kraken2': ['v_kraken2.txt', r"Kraken version (\S+)-beta"],
    'Quast': ['v_quast.txt', r"QUAST v(\S+)"],
    'CAT': ['v_cat.txt', r"CAT v(\S+)"]
}
results = OrderedDict()
results['nf-core/mag'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['fastp'] = '<span style="color:#999999;\">N/A</span>'
results['MEGAHIT'] = '<span style="color:#999999;\">N/A</span>'
results['Metabat'] = '<span style="color:#999999;\">N/A</span>'
results['NanoPlot'] = '<span style="color:#999999;\">N/A</span>'
results['Filtlong'] = '<span style="color:#999999;\">N/A</span>'
results['Porechop'] = '<span style="color:#999999;\">N/A</span>'
results['NanoLyse'] = '<span style="color:#999999;\">N/A</span>'
results['SPAdes'] = '<span style="color:#999999;\">N/A</span>'
results['BUSCO'] = '<span style="color:#999999;\">N/A</span>'
results['Centrifuge'] = '<span style="color:#999999;\">N/A</span>'
results['Kraken2'] = '<span style="color:#999999;\">N/A</span>'
results['CAT'] = '<span style="color:#999999;\">N/A</span>'
results['Quast'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/mag Software Versions'
section_href: 'https://github.com/nf-core/mag'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
