# nf-core/mag: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mag/usage](https://nf-co.re/mag/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

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

See the [nf-core/mag website documentation](https://nf-co.re/mag/usage#usage) for more information about pipeline specific parameters.

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

To allow also reproducible bin QC with BUSCO, the parameter `--save_busco_reference` can be used to store the reference database. This may be useful since BUSCO datasets are frequently updated and old versions do not always remain accessible.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Conda) - see below.

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
* `podman`
  * A generic configuration profile to be used with [Podman](https://podman.io/)
  * Pulls software from Docker Hub: [`nfcore/mag`](https://hub.docker.com/r/nfcore/mag/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity or Podman.
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

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

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

## Input specifications

The input data can be passed to nf-core/mag in two possible ways using the `--input` parameter.

### Direct FASTQ input (short reads only)

The easiest way is to specify directly the path (with wildcards) to your input FASTQ files. For example:

```bash
--input 'path/to/data/sample_*_R{1,2}.fastq.gz'
```

This input method only works with short read data and will assign all files to the same group. By default, this group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).

Please note the following additional requirements:

* Files names must be unique
* Valid file extensions: `.fastq.gz`, `.fastq`, `.fq.gz`, `.fq`
* The path must be enclosed in quotes
* The path must have at least one `*` wildcard character
* When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs
* To run single-end data you must additionally specify `--single_end`
* If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### Samplesheet input file

Alternatively, to assign different groups or to include long reads for hybrid assembly with metaSPAdes, you can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata.

This CSV file should contain the following columns:

`sample,group,short_reads_1,short_reads_2,long_reads`

The path to `long_reads` and `short_reads_2` is optional. Valid examples could look like the following:

```bash
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
```

or

```bash
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1.fastq.gz,,
sample2,0,data/sample2.fastq.gz,,
```

Please note the following requirements:

* 5 comma-seperated columns
* Valid file extension: `.csv`
* Must contain the header `sample,group,short_reads_1,short_reads_2,long_reads`
* Sample IDs must be unique
* `long_reads` can only be provided in combination with paired-end short read data
* Within one samplesheet either only single-end or only paired-end reads can be specified
* If single-end reads are specified, the command line parameter `--single_end` must be specified as well

Again, by default, the group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).
