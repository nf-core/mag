# nf-core/mag: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Quality control](#quality-control) of input reads - trimming and contaminant removal
* [Taxonomic classification](#taxonomic-classification) of trimmed reads
* [Assembly](#assembly) of trimmed reads
* [Binning](#binning) of assembled contigs
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Quality control

These steps trim away the adapter sequences present in input reads, trims away bad quality bases and sicard reads that are too short.
It also removes host contaminants and sequencing controls, such as PhiX or the Lambda phage.
FastQC is run for visualising the general quality metrics of the sequencing runs before and after trimming.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `QC_shortreads/fastqc`**

* `[sample]_R[1/2]_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files
* `[sample]_R[1/2].trimmed_fastqc.html`: FastQC report, containing quality metrics for trimmed and, if specified, filtered read files

### fastp

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the results folder and part of the MultiQC report.

**Output directory: `QC_shortreads/fastp/[sample]`**

* `fastp.html`: Interactive report
* `fastp.json`: Report in json format

### Remove PhiX sequences from short reads

The pipeline uses bowtie2 to map the reads against PhiX and removes mapped reads.

**Output directory: `QC_shortreads/remove_phix`**

* `[sample]_remove_phix_log.txt`: Contains a brief log file indicating how many reads have been retained.

### Host read removal

The pipeline uses bowtie2 to map short reads against the host reference genome specified with `--host_genome` or `--host_fasta` and removes mapped reads. The information about discarded and retained reads is also included in the MultiQC report.

**Output directory: `QC_shortreads/remove_host`**

* `[sample].bowtie2.log`: Contains the bowtie2 log file indicating how many reads have been mapped as well as a file listing the read ids of discarded reads.

### Remove Phage Lambda sequences from long reads

The pipeline uses Nanolyse to map the reads against the Lambda phage and removes mapped reads.

**Output directory: `QC_longreads/NanoLyse`**

* `[sample]_nanolyse_log.txt`: Contains a brief log file indicating how many reads have been retained.

### Filtlong and porechop

The pipeline uses filtlong and porechop to perform quality control of the long reads that are eventually provided with the `--manifest` option.

No direct host read removal is performed for long reads.
However, since within this pipeline filtlong uses a read quality based on k-mer matches to the already filtered short reads, reads not overlapping those short reads might be discarded.
The lower the parameter `--longreads_length_weight`, the higher the impact of the read qualities for filtering.
For further documentation see the [filtlong online documentation](https://github.com/rrwick/Filtlong).

### Quality visualisation for long reads

NanoPlot is used to calculate various metrics and plots about the quality and length distribution of long reads. For more information about NanoPlot see the [online documentation](https://github.com/wdecoster/NanoPlot).

**Output directory: `QC_longreads/NanoPlot_[sample]`**

* `raw_*.[png/html/txt]`: Plots and reports for raw data
* `filtered_*.[png/html/txt]`: Plots and reports for filtered data

## Taxonomic Classification

### Kraken

Kraken2 classifies reads using a k-mer based approach as well as assigns taxonomy using a Lowest Common Ancestor (LCA) algorithm.

**Output directory: `Taxonomy/kraken2/[sample]`**

* `kraken2.report`: Classification in the Kraken report format. See the [kraken manual](http://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports) for more details
* `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

### Centrifuge

Centrifuge is commonly used for the classification of DNA sequences from microbial samples. It uses an indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index.

More information on the [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) website

**Output directory: `Taxonomy/centrifuge/[sample]`**

* `report.txt`: Tab-delimited result file. See the [centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output) for information about the fields
* `kreport.txt`: Classification in the Kraken report format. See the [kraken manual](http://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports) for more details
* `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

### CAT

[CAT](https://github.com/dutilh/CAT) is a toolkit for annotating contigs and bins from metagenome-assembled-genomes. The MAG pipeline uses CAT to assign taxonomy to the contigs from megahit and/or SPAdes, and to assign taxonomy to genome bins based on the taxnomy of the contigs.

**Output directory: `Taxonomy/[assembler]`**

* `[assembler]-[sample].ORF2LCA.txt`: Tab-delimited files containing the lineage of each contig
* `[assembler]-[sample].names.txt`: Taxonomy classification, with names of each lineage levels instead og taxids
* `[assembler]-[sample].predicted_proteins.faa`: predicted protein sequences for each genome bins, in fasta format
* `[assembler]-[sample].predicted_proteins.gff`: predicted protein features for each genome bins, in gff format
* `[assembler]-[sample].log`: Log files
* `[assembler]-[sample].bin2classification.txt`: Taxonomy classification of the genome bins

## Assembly

Trimmed (short) reads are assembled with both megahit and SPAdes. Hybrid assembly is only supported by SPAdes.

### MEGAHIT

[MEGAHIT](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics short reads.

**Output directory: `Assembly/MEGAHIT`**

* `[sample].contigs.fa.gz`: Compressed metagenome assembly in fasta format
* `[sample].log`: Log file
* `[sample]_QC/`: Directory containing QUAST files

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) was originally a single genome assembler that later added support for assembling metagenomes.

**Output directory: `Assembly/SPAdes`**

* `[sample]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
* `[sample]_graph.gfa.gz`: Compressed assembly graph in gfa format
* `[sample]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
* `[sample]_log.txt`: Log file
* `[sample]_QC/`: Directory containing QUAST files

### SPAdesHybrid

SPAdesHybrid is a part of the [SPAdes](http://cab.spbu.ru/software/spades/) software and is used when the user provides both long and short reads.

**Output directory: `Assembly/SPAdesHybrid`**

* `[sample]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
* `[sample]_graph.gfa.gz`: Compressed assembly graph in gfa format
* `[sample]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
* `[sample]_log.txt`: Log file
* `[sample]_QC/`: Directory containing QUAST files

### Metagenome QC with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates metagenome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

**Output directory: `Assembly/[assembler]/[sample]_QC`**

* `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
* `quast.log`: QUAST log file
* `predicted_genes/[assembler]-[sample].rna.gff`: Contig positions for rRNA genes in gff version 3 format

## Binning

### Contig sequencing depth

Sequencing depth per contig and sample is generated by `jgi_summarize_bam_contig_depths --outputDepth`. The values correspond to `(sum of exactely aligned bases) / ((contig length)-2*75)`. For example, for two reads aligned exactly with `10` and `9` bases on a 1000 bp long contig the depth is calculated by `(10+9)/(1000-2*75)` (1000bp length of contig minus 75bp from each end, which is excluded).

**Output directory: `GenomeBinning`**

* `[assembler]-[sample]-depth.txt.gz`: Sequencing depth for each contig and sample, only for short reads.

### MetaBAT2

[MetaBAT2](https://bitbucket.org/berkeleylab/metabat) recovers genome bins (that is, contigs/scaffolds that all belongs to a same organism) from metagenome assemblies.

**Output directory: `GenomeBinning/MetaBAT2`**

* `[assembler]-[sample].*.fa`: Genome bins retrieved from input assembly
* `[assembler]-[sample].unbinned.*.fa`: Contigs that were not binned with other contigs but considered interesting. By default, these are at least 1 Mbp (`--min_length_unbinned_contigs`) in length and at most the 100 longest contigs (`--max_unbinned_contigs`) are reported

All the files and contigs in this folder will be assessed by QUAST and BUSCO.

**Output directory: `GenomeBinning/MetaBAT2/discarded`**

* `*.lowDepth.fa`: Low depth contigs that are filtered by MetaBat2
* `*.tooShort.fa`: Too short contigs that are filtered by MetaBat2
* `*.unbinned.pooled.fa`: Pooled unbinned contigs equal or above `--min_contig_size`, by default 1500 bp.
* `*.unbinned.remaining.fa`: Remaining unbinned contigs below `--min_contig_size`, by default 1500 bp, but not in any other file.

All the files in this folder contain small and/or unbinned contigs that are not further processed.

Files in these two folders contain all contigs of an assembly.

### QC for metagenome assembled genomes with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates genome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

**Output directory: `GenomeBinning/QC/QUAST/[assembler]-[bin]`**

* `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
* `quast.log`: QUAST log file
* `predicted_genes/[assembler]-[sample].rna.gff`: Contig positions for rRNA genes in gff version 3 format

**Output directory: `GenomeBinning/QC`**

* `quast_summary.tsv`: QUAST output for all bins summarized
* `quast_and_busco_summary.tsv`: Summary of BUSCO and QUAST results

### QC for metagenome assembled genomes with BUSCO

[BUSCO](https://busco-archive.ezlab.org/v3/) is a tool used to assess the completeness of a genome assembly. It is run on all the genome bins and high quality contigs obtained by MetaBAT2.

**Output directory: `GenomeBinning/QC/BUSCO`**

* `[assembler]-[bin]_busco_log.txt`: BUSCO log file
* `[assembler]-[bin]_busco.fna`: Nucleotide sequence of all identified BUSCOs
* `[assembler]-[bin]_busco.faa`: Aminoacid sequence of all identified BUSCOs

**Output directory: `GenomeBinning/QC`**

* `busco_summary.txt`: A summary table of the BUSCO results, with % of marker genes found
* `quast_and_busco_summary.tsv`; Summary of BUSCO and QUAST results

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output directory: `multiqc`**

* `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
* `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
* `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output directory: `pipeline_info`**

* Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
* Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
* Documentation for interpretation of results in HTML format: `results_description.html`.
