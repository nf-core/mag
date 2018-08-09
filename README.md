# [UNDER CONSTRUCTION] nf-core/mag

THIS PIPELINE IS A WORK IN PROGRESS. Thanks for checking it out! Hopefully it will be fonctional soon.

Assembly, binning and annotation of metagenomes

[![Build Status](https://travis-ci.org/nf-core/mag.svg?branch=master)](https://travis-ci.org/nf-core/mag)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/hadrieng/mag.svg)](https://hub.docker.com/r/hadrieng/mag)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
nf-core/mag: Assembly, binning and annotation of metagenomes

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The nf-core/mag pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline was written by [Hadrien Gourlé](https://hadriengourle.com) at [SLU](https://slu.se).
