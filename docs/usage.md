# nf-core/mag: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mag/usage](https://nf-co.re/mag/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Input specifications

The input data can be passed to nf-core/mag in two possible ways using the `--input` parameter.

### Direct FASTQ input (short reads only)

The easiest way is to specify directly the path (with wildcards) to your input FASTQ files. For example:

```console
--input 'path/to/data/sample_*_R{1,2}.fastq.gz'
```

This input method only works with short read data and will assign all files to the same group. By default, this group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).

Please note the following additional requirements:

- Files names must be unique
- Valid file extensions: `.fastq.gz`, `.fq.gz` (files must be compressed)
- The path must be enclosed in quotes
- The path must have at least one `*` wildcard character
- When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs
- To run single-end data you must additionally specify `--single_end`
- If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### Samplesheet input file

Alternatively, to assign different groups or to include long reads for hybrid assembly with metaSPAdes, you can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata.

This CSV file should contain the following columns:

`sample,group,short_reads_1,short_reads_2,long_reads`

The path to `long_reads` and `short_reads_2` is optional. Valid examples could look like the following:

```console
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
```

or

```console
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1.fastq.gz,,
sample2,0,data/sample2.fastq.gz,,
```

Please note the following requirements:

- 5 comma-seperated columns
- Valid file extension: `.csv`
- Must contain the header `sample,group,short_reads_1,short_reads_2,long_reads`
- Sample IDs must be unique
- FastQ files must be compressed (`.fastq.gz`, `.fq.gz`)
- `long_reads` can only be provided in combination with paired-end short read data
- Within one samplesheet either only single-end or only paired-end reads can be specified
- If single-end reads are specified, the command line parameter `--single_end` must be specified as well

Again, by default, the group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/mag --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTIDR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

See the [nf-core/mag website documentation](https://nf-co.re/mag/usage#usage) for more information about pipeline specific parameters.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/mag
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mag releases page](https://github.com/nf-core/mag/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

Additionally, to enable also reproducible results from the individual assembly tools this pipeline provides extra parameters. SPAdes is designed to be deterministic for a given number of threads. To generate reproducible results set the number of cpus with `--spades_fix_cpus` or `--spadeshybrid_fix_cpus`. This will overwrite the number of cpus specified in the `base.config` file and additionally ensure that it is not increased in case of retries for individual samples. MEGAHIT only generates reproducible results when run single-threaded.
You can fix this by using the prameter `--megahit_fix_cpu_1`. In both cases, do not specify the number of cpus for these processes in additional custom config files, this would result in an error.

MetaBAT2 is run by default with a fixed seed within this pipeline, thus producing reproducible results.

To allow also reproducible bin QC with BUSCO, run BUSCO providing already downloaded lineage datasets with `--busco_download_path` (BUSCO will be run using automated lineage selection in offline mode) or provide a specific lineage dataset via `--busco_reference` and use the parameter `--save_busco_reference`. This may be useful since BUSCO datasets are frequently updated and old versions do not always remain (easily) accessible.

For the taxonomic bin classification with [CAT](https://github.com/dutilh/CAT), when running the pipeline with `--cat_db_generate` the parameter `--save_cat_db` can be used to also save the generated database to allow reproducibility in future runs. Note that when specifying a pre-built database with `--cat_db`, currently the database can not be saved.

The taxonomic classification of bins with GTDB-Tk is not guaranteed to be reproducible, since the placement of bins in the reference tree is non-deterministic. However, the authors of the GTDB-Tk article examined the reproducibility on a set of 100 genomes across 50 trials and did not observe any difference (see [https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)).

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`, `test_hybrid`, `test_host_rm`, `test_hybrid_host_rm`, `test_busco_auto`
  - Profiles with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

Note, do not change number of CPUs with custom config files for the processes `spades`, `spadeshybrid` or `megahit` when specifying the parameters `--spades_fix_cpus`, `--spadeshybrid_fix_cpus` and `--megahit_fix_cpu_1` respectively.

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```

## A note on the ancient DNA subworkflow

MAG integrates an additional subworkflow to validate ancient DNA _de novo_ assembly:

[Characteristic patterns of ancient DNA (aDNA) damage](<(https://doi.org/10.1073/pnas.0704665104)>), namely DNA fragmentation and cytosine deamination (observed as C-to-T transitions) are typically used to authenticate aDNA sequences. By identifying assembled contigs carrying typical aDNA damages using [PyDamage](https://github.com/maxibor/pydamage), MAG can report and distinguish ancient contigs from contigs carrying no aDNA damage. Furthermore, to mitigate the effect of aDNA damage on contig sequence assembly, [freebayes](https://github.com/freebayes/freebayes) in combination with [BCFtools](https://github.com/samtools/bcftools) are used to (re)call the variants from the reads aligned to the contigs, and (re)generate contig consensus sequences.

## A note on bin refinement

DAS_Tool may not always be able to refine bins due to insufficient recovery of enough single-copy genes. In these cases you will get a NOTE such as

```console
[16/d330a6] NOTE: Process `NFCORE_MAG:MAG:BINNING_REFINEMENT:DASTOOL_DASTOOL (test_minigut_sample2)` terminated with an error exit status (1) -- Error is ignored
```

In this case, DAS_Tool pipeline has not necessarily failed but was unable to complete the refinement. You will therefore not expect to find any output files in the BinRefinement results directory for that sample.

If you are regularly getting such errors, you can try reducing the `--refine_bins_dastool_threshold` value, which will modify the scoring threshold defined in the [DAS_Tool publication](https://www.nature.com/articles/s41564-018-0171-1).

Furthermore, only
