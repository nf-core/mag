# ![nf-core/mag](docs/images/nf-core-mag_logo.png)

**Assembly, binning and annotation of metagenomes**.

[![GitHub Actions CI Status](https://github.com/nf-core/mag/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/mag/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/mag/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/mag/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/mag.svg)](https://hub.docker.com/r/nfcore/mag)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3589527.svg)](https://doi.org/10.5281/zenodo.3589527)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mag-4A154B?logo=slack)](https://nfcore.slack.com/channels/mag)

## Introduction

**nf-core/mag** is a bioinformatics best-practise analysis pipeline for assembly, binning, and annotation of metagenomes.

<p align="center">
    <img src="docs/images/mag_workflow.png" alt="nf-core/mag workflow overview" width="60%">
</p>

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/mag -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/mag -profile <docker/singularity/podman/conda/institute> --input '*_R{1,2}.fastq.gz'
    ```

    or

    ```bash
    nextflow run nf-core/mag -profile <docker/singularity/podman/conda/institute> --input 'samplesheet.csv'
    ```

See [usage docs](https://nf-co.re/mag/usage) and [parameter docs](https://nf-co.re/mag/parameters) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following: it supports both short and long reads, quality trims the reads and adapters with [fastp](https://github.com/OpenGene/fastp) and [Porechop](https://github.com/rrwick/Porechop), and performs basic QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
The pipeline then:

* assigns taxonomy to reads using [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) and/or [Kraken2](https://github.com/DerrickWood/kraken2/wiki)
* performs assembly using [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](http://cab.spbu.ru/software/spades/), and checks their quality using [Quast](http://quast.sourceforge.net/quast)
* performs metagenome binning using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/), and checks the quality of the genome bins using [Busco](https://busco.ezlab.org/)
* assigns taxonomy to bins using [CAT](https://github.com/dutilh/CAT)

Furthermore, the pipeline creates various reports in the results directory specified, including a [MultiQC](https://multiqc.info/) report summarizing some of the findings and software versions.

## Documentation

The nf-core/mag pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/mag/usage) and [output](https://nf-co.re/mag/output). Detailed information about how to specify the input can be found under [input specifications](https://nf-co.re/mag/usage#input_specifications).

### Group-wise co-assembly and co-abundance computation

Each sample has an associated group ID (see [input specifications](https://nf-co.re/mag/usage#input_specifications)). This group information can be used for group-wise co-assembly with `MEGAHIT` or `SPAdes` and/or to compute co-abundances for the binning step with `MetaBAT2`. By default, group-wise co-assembly is disabled, while the computation of group-wise co-abundances is enabled. For more information about how this group information can be used see the documentation for the parameters [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode).

When group-wise co-assembly is enabled, `SPAdes` is run on accordingly pooled read files, since `metaSPAdes` does not yet allow the input of multiple samples or libraries. In contrast, `MEGAHIT` is run for each group while supplying lists of the individual readfiles.

## Credits

nf-core/mag was written by [Hadrien Gourlé](https://hadriengourle.com) at [SLU](https://slu.se), [Daniel Straub](https://github.com/d4straub) and [Sabrina Krakau](https://github.com/skrakau) at the [Quantitative Biology Center (QBiC)](http://qbic.life).

Long read processing was inspired by [caspargross/HybridAssembly](https://github.com/caspargross/HybridAssembly) written by Caspar Gross [@caspargross](https://github.com/caspargross)

Many thanks to the additional contributors who have helped out and/or provided suggestions:

* [Alexander Peltzer](https://github.com/apeltzer)
* [Phil Ewels](https://github.com/ewels)
* [Gisela Gabernet](https://github.com/ggabernet)
* [Harshil Patel](https://github.com/drpatelh)
* [Johannes Alneberg](https://github.com/alneberg)
* [Maxime Borry](https://github.com/maxibor)
* [Maxime Garcia](https://github.com/MaxUlysse)
* [Michael L Heuer](https://github.com/heuermh)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#mag` channel](https://nfcore.slack.com/channels/mag) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/mag for your analysis, please cite it using the following doi: [10.5281/zenodo.3589527](https://doi.org/10.5281/zenodo.3589527)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

References of tools used in this pipeline are as follows:

* **Bowtie2**  Langmead, B. and Salzberg, S. L. 2012 Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), p. 357–359. doi: [10.1038/nmeth.1923](https:/dx.doi.org/10.1038/nmeth.1923).
* **Busco** Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUSCO: assessing genome assembly and annotation completeness. In Gene prediction (pp. 227-245). Humana, New York, NY. [https://doi.org/10.1007/978-1-4939-9173-0_14](https://doi.org/10.1007/978-1-4939-9173-0_14)
* **CAT** von Meijenfeldt, F. B., Arkhipova, K., Cambuy, D. D., Coutinho, F. H., & Dutilh, B. E. (2019). Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT. Genome biology, 20(1), 1-14. [https://doi.org/10.1186/s13059-019-1817-x](https://doi.org/10.1186/s13059-019-1817-x). Home: [https://github.com/dutilh/CAT](https://github.com/dutilh/CAT)
* **Centrifuge** Kim, D., Song, L., Breitwieser, F. P., & Salzberg, S. L. (2016). Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome research, 26(12), 1721-1729. [https://doi.org/10.1101/gr.210641.116](https://doi.org/10.1101/gr.210641.116). Home: [http://www.ccb.jhu.edu/software/centrifuge](http://www.ccb.jhu.edu/software/centrifuge)
* **FastP** Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics , 34(17), i884–i890. [https://doi.org/10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560). Home: [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)
* **FastQC** Home: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* **Filtlong** Home: [https://github.com/rrwick/Filtlong](https://github.com/rrwick/Filtlong)
* **Kraken2** Wood, D et al., 2019. Improved metagenomic analysis with Kraken 2. Genome Biology volume 20, Article number: 257. [https://doi.org/10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0). Home: [https://ccb.jhu.edu/software/kraken2/](https://ccb.jhu.edu/software/kraken2/)
* **Krona** Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC bioinformatics, 12(1), 1-10. [https://doi.org/10.1186/1471-2105-12-385](https://doi.org/10.1186/1471-2105-12-385). Home: [https://github.com/marbl/Krona/wiki](https://github.com/marbl/Krona/wiki).
* **MEGAHIT** Li, D., Luo, R., Liu, C. M., Leung, C. M., Ting, H. F., Sadakane, K., ... & Lam, T. W. (2016). MEGAHIT v1. 0: a fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods, 102, 3-11. [https://doi.org/10.1016/j.ymeth.2016.02.020](https://doi.org/10.1016/j.ymeth.2016.02.020). Home: [https://github.com/voutcn/megahit](https://github.com/voutcn/megahit)
* **MetaBAT2** Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359. [https://doi.org/10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359). Home: [https://bitbucket.org/berkeleylab/metabat](https://bitbucket.org/berkeleylab/metabat)
* **MultiQC** Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354). Home: [https://multiqc.info/](https://multiqc.info/)
* **NanoLyse** De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666-2669. [https://doi.org/10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149). Home: [https://github.com/wdecoster/nanolyse](https://github.com/wdecoster/nanolyse)
* **NanoPlot** De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666-2669. [https://doi.org/10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149). Home: [https://github.com/wdecoster/NanoPlot](https://github.com/wdecoster/NanoPlot)
* **Porechop** Home: [https://github.com/rrwick/Porechop](https://github.com/rrwick/Porechop)
* **SAMtools** Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics , 25(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352). Home: [http://www.htslib.org/](http://www.htslib.org/)
* **SPAdes** Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, 27(5), 824-834. [https://doi.org/10.1101/gr.213959.116](https://doi.org/10.1101/gr.213959.116)

In addition, this repository uses test data from the following study:

* Bertrand, D., Shaw, J., Kalathiyappan, M., Ng, A. H. Q., Kumar, M. S., Li, C., ... & Nagarajan, N. (2019). Hybrid metagenomic assembly enables high-resolution analysis of resistance determinants and mobile elements in human microbiomes. Nature biotechnology, 37(8), 937-944. [https://doi.org/10.1038/s41587-019-0191-2](https://doi.org/10.1038/s41587-019-0191-2)
