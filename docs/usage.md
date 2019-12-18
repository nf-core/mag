# nf-core/mag: Usage

## Table of contents

- [Introduction](#general-nextflow-info)
- [Running the pipeline](#running-the-pipeline)
- [Updating the pipeline](#updating-the-pipeline)
- [Reproducibility](#reproducibility)
- [Main arguments](#main-arguments)
  - [`-profile`](#-profile-single-dash)
    - [`docker`](#docker)
    - [`awsbatch`](#awsbatch)
    - [`standard`](#standard)
    - [`none`](#none)
  - [`--reads`](#--reads)
  - [`--singleEnd`](#--singleend)
  - [`--manifest`](#--manifest)
- [Optional Arguments](#optional-arguments)
  - [Trimming Options](#trimming-options)
  - [Trimming Options for long reads](#trimming-options-for-long-reads)
  - [Taxonomic classification](#taxonomic-classification)
  - [Binning Options](#binning-options)
- [Job Resources](#job-resources)
- [Automatic resubmission](#automatic-resubmission)
- [Custom resource requests](#custom-resource-requests)
- [AWS batch specific parameters](#aws-batch-specific-parameters)
  - [`-awsbatch`](#-awsbatch)
  - [`--awsqueue`](#--awsqueue)
  - [`--awsregion`](#--awsregion)
- [Other command line parameters](#other-command-line-parameters)
  - [`--outdir`](#--outdir)
  - [`--email`](#--email)
  - [`-name`](#-name-single-dash)
  - [`-resume`](#-resume-single-dash)
  - [`-c`](#-c-single-dash)
  - [`--max_memory`](#--max_memory)
  - [`--max_time`](#--max_time)
  - [`--max_cpus`](#--max_cpus)
  - [`--plaintext_emails`](#--plaintext_emails)
  - [`--sampleLevel`](#--sampleLevel)
  - [`--multiqc_config`](#--multiqc_config)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/mag --reads '*_R{1,2}.fastq.gz' -profile docker
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

First, go to the [nf-core/mag releases page](https://github.com/nf-core/mag/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

- `awsbatch`
  - A generic configuration profile to be used with AWS Batch.
- `conda`
  - A generic configuration profile to be used with [conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from dockerhub: [`nfcore/mag`](http://hub.docker.com/r/nfcore/mag/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  - Pulls software from DockerHub: [`nfcore/mag`](http://hub.docker.com/r/nfcore/mag/)
- `test`, `test_hybrid`
  - Profiles with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--manifest`

The pipeline has support for hybrid (with long and short reads) assembly, with the `--manifest` option.
The option take a tab-separated file with 4 headerless columns: Sample_Id, Long_Reads, Short_Reads_1, Short_Reads_2
Only one file path per entry is allowed, and single-end short reads are not supported.

## Trimming options

### `--adapter_forward`

Sequence of 3' adapter to remove in the forward reads

### `--adapter_reverse`

Sequence of 3' adapter to remove in the reverse reads

### `--mean_quality`

Mean qualified quality value for keeping read (default: 15)

### `--trimming_quality`

Trimming quality value for the sliding window (default: 15)

### `--keep_phix`

Keep reads similar to the Illumina internal standard PhiX genome (default: false)

## Trimming options for long reads

### `--skip_adapter_trimming`

Skip removing adapter sequences from long reads

### `--longreads_min_length`

Discard any read which is shorter than this value (default: 1000)

### `--longreads_keep_percent`

Keep this percent of bases (default: 90)

### `--longreads_length_weight`

The higher the more important is read length when choosing the best reads (default: 10)

### `--keep_lambda`

Keep reads similar to the ONT internal standard Escherichia virus Lambda genome (default: false)

## Taxonomic classification

Taxonomic classification is disabled by default.
You have to specify one of the options below to activate it.

### `--centrifuge_db`

Database for taxonomic binning with centrifuge (default: none). E.g. "<ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz>"

### `--kraken2_db`

Database for taxonomic binning with kraken2 (default: none). E.g. "<ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz>"

### `--cat_db`

Database for taxonomic classification of metagenome assembled genomes (default: none). E.g. "<tbb.bio.uu.nl/bastiaan/CAT*prepare/CAT_prepare_20190108.tar.gz>"
The zipped file needs to contain a folder named "\_taxonomy*" and "_CAT_database_" that hold the respective files.

## Binning options

### `--min_contig_size`

Minimum contig size to be considered for binning (default: 1500)

### `--busco_reference`

Download path for BUSCO database, available databases are listed here: <https://busco.ezlab.org/>
(default: <https://busco.ezlab.org/datasets/bacteria_odb9.tar.gz>)

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-core-invite.herokuapp.com/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
