#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/mag': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'fastqc': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'fastp': ['v_fastp.txt', r"fastp (\S+)"],
    'megahit': ['v_megahit.txt', r"MEGAHIT v(\S+)"]
}
results = OrderedDict()
results['nf-core/mag'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['fastqc'] = '<span style="color:#999999;\">N/A</span>'
results['fastp'] = '<span style="color:#999999;\">N/A</span>'
results['megahit'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print('''
id: 'nf-core/mag-software-versions'
section_name: 'nf-core/mag Software Versions'
section_href: 'https://github.com/nf-core/mag'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k, v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k, v))
print("    </dl>")
