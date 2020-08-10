# nf-core/mag: Usage

## Table of contents

* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Core Nextflow arguments](#core-nextflow-arguments)
  * [`-profile`](#-profile)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
    * [`Custom resource requests`](#custom-resource-requests)
  * [`Running in the background`](#running-in-the-background)
    * [`Nextflow memory requirements`](#nextflow-memory-requirements)
* [Pipeline specific arguments](#pipeline-specific-arguments)
  * [`--input`](#--input)
  * [`--single_end`](#--single_end)
  * [`--manifest`](#--manifest)
* [Quality control (for short reads)](#quality-control-for-short-reads)
  * [`--adapter_forward`](#--adapter_forward)
  * [`--adapter_reverse`](#--adapter_reverse)
  * [`--mean_quality`](#--mean_quality)
  * [`--trimming_quality`](#--trimming_quality)
  * [`--host_fasta`](#--host_fasta)
  * [`--host_genome (using iGenomes)`](#--host_genome-using-iGenomes)
  * [`--host_removal_verysensitive`](#--host_removal_verysensitive)
  * [`--keep_phix`](#--keep_phix)
* [Quality control for long reads](#quality-control-for-long-reads)
  * [`--skip_adapter_trimming`](#--skip_adapter_trimming)
  * [`--longreads_min_length`](#--longreads_min_length)
  * [`--longreads_keep_percent`](#--longreads_keep_percent)
  * [`--longreads_length_weight`](#--longreads_length_weight)
  * [`--keep_lambda`](#--keep_lambda)
* [Taxonomic classification](#taxonomic-classification)
  * [`--centrifuge_db`](#--centrifuge_db)
  * [`--kraken2_db`](#--kraken2_db)
  * [`--cat_db`](#--cat_db)
* [Assembly options](#assembly-options)
  * [`--megahit_fix_cpu_1`](#--megahit_fix_cpu_1)
  * [`--spades_fix_cpus`](#--spades_fix_cpus)
  * [`--spadeshybrid_fix_cpus`](#--spadeshybrid_fix_cpus)
* [Binning Options](#binning-options)
  * [`--min_contig_size`](#--min_contig_size)
  * [`--busco_reference`](#--busco_reference)
  * [`--min_length_unbinned_contigs`](#--min_length_unbinned_contigs)
  * [`--max_unbinned_contigs`](#--max_unbinned_contigs)
  * [`--metabat_rng_seed`](#--metabat_rng_seed)

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/mag --input '*_R{1,2}.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/mag
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mag releases page](https://github.com/nf-core/mag/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

Additionally, to enable also reproducible results from the individual assembly tools this pipeline provides extra parameters. SPAdes is designed to be deterministic for a given number of threads. To generate reproducible results set the number of cpus with `--spades_fix_cpus` or `--spadeshybrid_fix_cpus`. This will overwrite the number of cpus specified in the `base.config` file and additionally ensure that it is not increased in case of retries for individual samples. MEGAHIT only generates reproducible results when run single-threaded.
You can fix this by using the prameter `--megahit_fix_cpu_1`. In both cases, do not specify the number of cpus for these processes in additional custom config files, this would result in an error.

MetaBAT2 is run by default with a fixed seed within this pipeline, thus producing reproducible results.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/mag`](https://hub.docker.com/r/nfcore/mag/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/mag`](https://hub.docker.com/r/nfcore/mag/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`, `test_hybrid`, `test_host_rm`, `test_hybrid_host_rm`
  * Profiles with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `spades` 100GB of memory, you could use the following config:

```nextflow
process {
  withName: spades {
    memory = 100.GB
  }
}
```

Do not change number of cpus for the processes `spades`, `spadeshybrid` or `megahit` in combination with the parameters `--spades_fix_cpus`, `--spadeshybrid_fix_cpus` and `--megahit_fix_cpu_1` respectively.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Pipeline specific arguments

### `--input`

Use this to specify the location of your input FastQ files. For example:

```bash
--input 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--input`. For example:

```bash
--single_end --input '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--manifest`

The pipeline has support for hybrid (with long and short reads) assembly, with the `--manifest` option.
The option take a tab-separated file with 4 headerless columns: Sample_Id, Long_Reads, Short_Reads_1, Short_Reads_2
Only one file path per entry is allowed, and single-end short reads are not supported.

## Quality control (for short reads)

### `--adapter_forward`

Sequence of 3' adapter to remove in the forward reads

### `--adapter_reverse`

Sequence of 3' adapter to remove in the reverse reads

### `--mean_quality`

Mean qualified quality value for keeping read (default: 15)

### `--trimming_quality`

Trimming quality value for the sliding window (default: 15)

### `--host_fasta`

Use this to speficify the full path to a host reference FASTA file to remove host reads with Bowtie 2 (default: none).

### `--host_genome (using iGenomes)`

Alternatively, you can specify an iGenomes reference that should be used for host read removal with Bowtie 2.

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

There are 31 different species supported in the iGenomes references. You can specify which genome to use with the `--host_genome` flag.
You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config).

When using the `--host_genome` parameter, both the iGenomes FASTA file as well as corresponding, already pre-built Bowtie 2 index files will be used for host read removal.

### `--host_removal_verysensitive`

Use the `--very-sensitive` setting (instead of `--sensitive`) for Bowtie 2 to map reads against host genome (default: false)

### `--keep_phix`

Keep reads similar to the Illumina internal standard PhiX genome (default: false)

## Quality control for long reads

### `--skip_adapter_trimming`

Skip removing adapter sequences from long reads

### `--longreads_min_length`

Discard any read which is shorter than this value (default: 1000)

### `--longreads_keep_percent`

Keep this percent of bases (default: 90)

### `--longreads_length_weight`

The higher the more important is read length when choosing the best reads (default: 10).
The default value focuses on length instead of quality to improve assembly size.
In order to assign equal weights to read lengths and read qualities set this parameter to 1.
This might be useful, for example, to benefit indirectly from the removal of short host reads (causing lower qualities for reads not overlapping filtered short reads).

### `--keep_lambda`

Keep reads similar to the ONT internal standard Escherichia virus Lambda genome (default: false)

## Taxonomic classification

Taxonomic classification is disabled by default.
You have to specify one of the options below to activate it.

### `--centrifuge_db`

Database for taxonomic binning with centrifuge (default: none). E.g. "<ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz>"

### `--kraken2_db`

Database for taxonomic binning with kraken2 (default: none). E.g. "<ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz>"

### `--cat_db`

Database for taxonomic classification of metagenome assembled genomes (default: none). E.g. "<http://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz>"
The zipped file needs to contain a folder named "\_taxonomy*" and "_CAT_database_" that hold the respective files.

## Assembly options

### `--megahit_fix_cpu_1`

Fix number of CPUs for MEGAHIT to 1. Not increased with retries (default: false). See also [Reproducibility](#reproducibility).

### `--spades_fix_cpus`

Fixed number of CPUs used by SPAdes. Not increased with retries (default: none).

### `--spadeshybrid_fix_cpus`

Fixed number of CPUs used by SPAdes hybrid. Not increased with retries (default: none).

## Binning options

### `--min_contig_size`

Minimum contig size to be considered for binning, for forwarding into downstream analysis (i.e. QUAST and BUSCO) and reporting (default: 1500)

### `--busco_reference`

Download path for BUSCO database, available databases are listed here: <https://busco.ezlab.org/>
(default: <https://busco.ezlab.org/datasets/bacteria_odb9.tar.gz>)

### `--min_length_unbinned_contigs`

Minimal length of contigs that are not part of any bin but treated as individual genome (default: 1000000)
Contigs that do not fulfill the thresholds of `--min_length_unbinned_contigs` and `--max_unbinned_contigs` are pooled for downstream analysis and reporting, except contigs that also do not fullfill `--min_contig_size` are not considered further.

### `--max_unbinned_contigs`

Maximal number of contigs that are not part of any bin but treated as individual genome (default: 100)
Contigs that do not fulfill the thresholds of `--min_length_unbinned_contigs` and `--max_unbinned_contigs` are pooled for downstream analysis and reporting, except contigs that also do not fullfill `--min_contig_size` are not considered further.

### `--metabat_rng_seed`

RNG seed for MetaBAT2. Use postive integer to ensure reproducibility (default: 1).
Set to 0 to use random seed.
