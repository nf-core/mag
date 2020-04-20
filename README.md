# ![mag](https://raw.githubusercontent.com/nf-core/mag/master/docs/images/mag_logo.png)

**Assembly, binning and annotation of metagenomes**.

[![GitHub Actions](https://github.com/nf-core/mag/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/mag/workflows/nf-core%20CI/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/mag.svg)](https://hub.docker.com/r/nfcore/mag)

## Introduction

This pipeline is for assembly, binning, and annotation of metagenomes.
It supports both short and long reads, quality trims the reads and adapters with [fastp](https://github.com/OpenGene/fastp) and [porechop](https://github.com/rrwick/Porechop), and performs basic QC with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

The pipeline then:

- assigns taxonomy to reads using [centrifuge](https://ccb.jhu.edu/software/centrifuge/) and/or [kraken2](https://ccb.jhu.edu/software/kraken2/)
- performs assembly using [megahit](https://github.com/voutcn/megahit) and [spades](http://cab.spbu.ru/software/spades/), and checks their quality using [quast](http://quast.sourceforge.net/quast)
- performs metagenome binning using [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/), and checks the quality of the genome bins using [busco](https://busco.ezlab.org/)

Furthermore, the pipeline creates various reports in the results directory specified, including a [multiqc](https://multiqc.info/) report summarizing some of the findings and software versions.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool that runs tasks across multiple compute infrastructures in a very portable manner. It comes with Docker containers making installation trivial and results highly reproducible.

## Documentation

The nf-core/mag pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
   - [Local installation](https://nf-co.re/usage/local_installation)
   - [Adding your own system config](https://nf-co.re/usage/adding_own_config)
   - [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

This pipeline was written by [Hadrien Gourlé](https://hadriengourle.com) at [SLU](https://slu.se) and Daniel Straub ([@d4straub](https://github.com/d4straub)).

Long read processing was inspired by [caspargross/HybridAssembly](https://github.com/caspargross/HybridAssembly) written by Caspar Gross [@caspargross](https://github.com/caspargross)
