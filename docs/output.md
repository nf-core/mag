# nf-core/mag: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

- [FastQC](#fastqc) - read quality control
- [fastp](#fastp) - read trimming
- [megahit](#megahit) - assembly
- [quast](#quast) - assembly quality report
- [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
- [metabat](#metabat) - recovered draft genome bins from assembled metagenomes
- [checkm](#checkm) - genome bins quality control and merging of compatible bins
- [refinem](#refinem) - improving metagenome-assembled genome bins

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows both _untrimmed_ and _trimmed_ reads.

**Output directory: `results/fastqc`**

- `sample_fastqc.html`
  - FastQC report, containing quality metrics for your untrimmed raw fastq files
- `zips/sample_fastqc.zip`
  - zip file containing the FastQC report, tab-delimited data file and plot images

## fastp

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads.

**Output directory: `None`**

- The trimmed reads are not included in the output

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

- `Project_multiqc_report.html`
  - MultiQC report - a standalone HTML file that can be viewed in your web browser
- `Project_multiqc_data/`
  - Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info

## Megahit

[megahit](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics HTS reads.

**Output directory: `results/megahit`**

- `sample.fasta`
  - metagenome assembly in fasta format

## Quast

[quast](http://cab.spbu.ru/software/quast/) is a tool that evaluates genome and metagenome assemblies by computing various metrics

**output directory:** `None`

- The quast output is included in the multiqc report

## Metabat

[metabat](https://bitbucket.org/berkeleylab/metabat) recovers genome bins (that is, contigs/scaffolds that all bolongs to a same organism) from metagenome assemblies.

**output directory: `results/metabat`**

- `sample.bam`
  - reads mapped against the megahit assembly
- `bins/sample_X.fa`
  - the putative genome bins retrieved by metabat

## checkm

[checkm](https://github.com/Ecogenomics/CheckM) Assess the quality of microbial genomes recovered from isolates, single cells, and metagenomes.

**output directory: `results/checkm`**

- `sample/`
  - directory containing the improved bins
- `sample_stats`
  - directory containing stats about completeness and contamination of the bins, as well as plots.

## refinem

_(this process is optional)_

[refinem](https://github.com/dparks1134/RefineM) is a companion tool to checkm, that refines and filter incongruent contigs from the genome bins.

**output directory: `results/refinem`**

- `sample.X.fa`
  - the genome bins, improved by refinem
