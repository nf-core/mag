# nf-core/mag: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mag/usage](https://nf-co.re/mag/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Input specifications

The input data can be passed to nf-core/mag in two possible ways, either using the `--input` parameter of raw-reads alone or `--input` additionally with `--assembly_input` that specifies pre-built assemblies.

### Samplesheet input file

You can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata. Furthermore when a `run` column is present, the pipeline will also perform run- or lane-wise concatenation, for cases where you may have a sample or library sequenced with the same sequencing configuration across multiple runs. The optional run merging happens after short read QC (adapter clipping, host/PhiX removal etc.), and prior to normalisation, taxonomic profiling, and assembly.

If short reads are provided (`short_reads_1`), the `short_reads_platform` column is required. Valid options include:

- `ILLUMINA`, `BGISEQ`, `LS454`, `ION_TORRENT`, `DNBSEQ`, `ELEMENT`, `ULTIMA`, `VELA_DIAGNOSTICS`, `GENAPSYS`, `GENEMIND`, `TAPESTRI`.

If long reads are provided (`long_reads`), the `long_reads_platform` column is required. Valid options include:

- `OXFORD_NANOPORE`: Oxford Nanopore Technologies (ONT) reads, which may have higher error rates compared to newer technologies.
- `OXFORD_NANOPORE_HQ`: High-quality ONT reads, typically with an error rate of less than 5%, achievable with the latest ONT chemistries and sequencing platforms. This option should generally be used unless working with older data.
- `PACBIO_SMRT`: Pacific Biosciences Single Molecule Real-Time (SMRT) sequencing reads.

These platform fields are important for downstream alignment and assembly tools.

An nf-core/mag input samplesheet file can contain the following columns:

`sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform,long_reads_platform`

In cases of short-read only data run, the path to `long_reads` and `short_reads_2` is optional. Valid examples could look like the following:

```csv title="samplesheet_mix.csv"
sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform,long_reads_platform
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz,ILLUMINA,OXFORD_NANOPORE
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz,ILLUMINA,OXFORD_NANOPORE
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,,ILLUMINA,
```

or

```csv title="samplesheet_shortreadonly.csv"
sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform
sample1,0,data/sample1.fastq.gz,,,ILLUMINA
sample2,0,data/sample2.fastq.gz,,,ILLUMINA
```

or to additionally perform run merging of two runs from sample1:

```csv title="samplesheet_mix_mergeruns.csv"
sample,run,group,short_reads_1,short_reads_2,long_reads,short_reads_platform,long_reads_platform
sample1,1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz,ILLUMINA,OXFORD_NANOPORE
sample1,2,0,data/sample1_lane2_R1.fastq.gz,data/sample1_lane2_R2.fastq.gz,data/sample1.fastq.gz,ILLUMINA,OXFORD_NANOPORE
sample2,0,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz,ILLUMINA,OXFORD_NANOPORE
sample3,1,0,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,,ILLUMINA,OXFORD_NANOPORE
```

If only long read data is available, the columns `short_reads_1` and `short_reads_2` can be left out:

```csv title="samplesheet_longreadonly.csv"
sample,run,group,long_reads,long_reads_platform
sample1,1,0,data/sample1A.fastq.gz,OXFORD_NANOPORE
sample1,2,0,data/sample1B.fastq.gz,OXFORD_NANOPORE
sample2,0,0,data/sample2.fastq.gz,OXFORD_NANOPORE
sample3,1,0,data/sample3.fastq.gz,OXFORD_NANOPORE
```

In this case only long-read only assemblies will be able to be executed (e.g. Flye or MetaMDBG).

Please note the following requirements:

- a minimum of 5 comma-separated columns
- Valid file extension: `.csv`
- Must contain the header `sample,group,short_reads_1,short_reads_2,long_reads` (where `run` can be optionally added)
- Run IDs must be unique within a multi-run sample. A sample with multiple runs will be automatically concatenated.
- FastQ files must be compressed (`.fastq.gz`, `.fq.gz`)
- Within one samplesheet either only single-end or only paired-end short reads can be specified
- If single-end reads are specified, the command line parameter `--single_end` must be specified as well

### The `group` column

In metagenomics, there are different strategies to deal with the _assembly_ and _binning_ of multiple samples.
These are usually referred to as "individual" or "single" assembly and binning, versus "co-" or "pooled" assembly and binning.
In the former, samples are processed independently.
In the latter, samples are processed together.
A common strategy is to assemble each sample individually (single assembly) and then pool all of the contigs together and map the reads against them (co-binning).
This is usually chosen since assembly is a computationally intensive process, it is very costly to assemble all samples together.
For binning, however, resources aren't as limiting, and binning algorithms can leverage the fact that there are multiple samples from which to draw information, which can improve the quality of output bins.

nf-core/mag, by default, follows this approach: the group information from the input sample sheet in is only used to compute co-abundances for the binning step (co-binning), but not for group-wise co-assembly (thus single assembly).
That means that if you define one group for all of your samples, they will be assembled individually, and then binned in a pooled fashion, with samples being mapped to all contigs of all other samples.

If you'd like to also _assemble_ your samples in a pooled fashion (co-assembly), see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode).

### Supplying pre-computed assemblies

It is also possible to run nf-core/mag on pre-computed assemblies, by supplying a CSV file to the parameter `--assembly_input` in addition to the raw reads supplied to `--input`. Supplying assembly input skips all read pre-processing and assembly, jumping straight to the binning stage of the pipeline.

The assembly CSV file should contain the following columns:

`id,group,assembler,fasta`

Where `id` is the ID of the assembly, group is the assembly/binning group (see samplesheet information section for more details), `assembler` is the assembler used to produce the assembly (one of `MEGAHIT`, `SPAdes`, `SPAdesHybrid`, `Flye` or `MetaMDBG`), and `fasta` is the path to the assembly fasta file. Input fasta files can be compressed or uncompressed, but compressed assemblies will be automatically uncompressed for use within the pipeline. The exact information required for each supplied assembly depends on whether the assemblies provided are single assemblies or group-wise co-assemblies. For the following example `--input` CSV:

```csv title="samplesheet.csv"
sample,group,short_reads_1,short_reads_2,short_reads_platform
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,ILLUMINA
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,ILLUMINA
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,ILLUMINA
```

If the assemblies are single assemblies, then the `id` and `group` columns should match those supplied in the `-input` read CSV files for each read set:

```csv title="samplesheet.csv"
id,group,assembler,fasta
sample1,0,MEGAHIT,MEGAHIT-sample1.contigs.fa.gz
sample1,0,SPAdes,SPAdes-sample1.fasta.gz
sample2,0,MEGAHIT,MEGAHIT-sample2.contigs.fa.gz
sample2,0,SPAdes,SPAdes-sample2.contigs.fasta.gz
sample3,1,MEGAHIT,MEGAHIT-sample3.contigs.fa.gz
sample3,1,SPAdes,SPAdes-sample3.contigs.fasta.gz
```

If the assemblies are co-assemblies, the parameter `--coassemble_group` should additionally be specified. In this case, the `id` column should uniquely identify the assembly, while `group` should match those specified in the `--input` CSV file:

```csv title="samplesheet.csv"
id,group,assembler,fasta
group-0,0,MEGAHIT,MEGAHIT-group-0.contigs.fa.gz
group-0,0,SPAdes,SPAdes-group-0.contigs.fasta.gz
group-1,1,MEGAHIT,MEGAHIT-group-1.contigs.fa.gz
group-1,1,SPAdes,SPAdes-group-1.contigs.fasta.gz
```

When supplying pre-computed assemblies, reads **must** also be provided in the CSV input format to `--input`, and should be the reads used to build the assemblies, i.e., adapter-removed, run-merged etc.. Preprocessing steps will not be ran on raw reads when pre-computed assemblies are supplied. As long reads are only used for assembly, any long read fastq files listed in the reads CSV are ignored.

### Databases

nf-core/mag uses multiple tools that require additional databases.
The pipeline will download these databases for you if not supplied via a pipeline-level parameter.
However, we generally recommend you download these databases manually once, and place them in a long-term storage directory or cache from where you can re-use the databases across runs.
Here we provide specific additional information for some of the downloading databases which can be more tricky to prepare.

For any other database, or in doubt, please check the [parameters page](https://nf-co.re/mag/parameters).

#### BUSCO

BUSCO can download lineage datasets automatically as it needs them, but this process can be slow if the internet connection is unstable, and may lead to redundant downloads across different samples.
To avoid this, you can provide a local copy of a BUSCO lineage dataset using the `--busco_db` parameter.
The local directory must follow a specific structure for BUSCO to recognize it.

To prepare the lineage dataset, the recommended method is to use BUSCO's built-in download functionality before running the pipeline:

```bash
busco --download_path <your_db> --download <lineage>
```

This command downloads the specified lineage into the `<your_db>/` directory and creates the required directory structure.
`<lineage>` can be any [supported dataset](https://busco-data.ezlab.org/v5/data/lineages/), such as `alphaproteobacteria_odb12`.
You can also specify `prokaryota` or `all` to download multiple lineages, which is necessary for automatic lineage selection in offline mode.

Alternatively, you can manually download a specific lineage tarball from [https://busco-data.ezlab.org/v5/data/lineages/](https://busco-data.ezlab.org/v5/data/lineages/), extract it, and place it in the appropriate location: `<your_db>/lineages/<taxa>_odb<XX>`.
The lineage directory (e.g., `bacteria_odb12`) should contain files such as `dataset.cfg` and a `hmms/` subdirectory at the top level.
Then, you must provide `--busco_db <your_db>/` and `--busco_db_lineage <downloaded_lineage>` to the pipeline.
You can also pass to `--busco_db` a URL pointing to a lineage tarball, or the tarball itself if stored locally.

> [!WARNING]
> When any kind of parameter is provided to `--busco_db`, BUSCO will run in offline mode.
> If the lineage specified via `--busco_db_lineage` is not found locally, or if you attempt automatic lineage selection without having a complete lineage dataset pre-downloaded, BUSCO will fail.

### CAT

> [!WARNING]
> The default database (CAT_nr) is very large at ~200 GB!
> This can take a long time, so we strongly recommend downloading and unzipping prior the pipeline run.

This database can be downloaded from the [CAT developers' website](https://tbb.bio.uu.nl/tina/CAT_pack_prepare/), which is based in the Netherlands (and could be slow for other regions of the world).

> [!NOTE]
> By default, the pipeline will use the `NCBI nr` database.

Enabling the `--cat_db_generate` option will create a new database using the latest version of the NCBI nr database.
This requires a large download (over 200 GB) and several hours of subsequent processing.
If you enable the `--save_cat_db` option, the database will be saved in the `Taxonomy/CAT/db` directory and can be reused in future runs by specifying its path with `--cat_db`.

### GTDB

> [!WARNING]
> This database is very large at ~110 GB!
> This can take a long time, so we strongly recommend downloading and unzipping prior the pipeline run.

This database can be downloaded from the [GTDB developers' website](default: https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/), which is based in Australia (and could be slow for other regions of the world).
The developers also offer a 'split' archive of 10GB files that can be downloaded more stably from [here](https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/split_package/) and subsequently (manually) combined after.
More documentation can be seen [here](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data).

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/mag --input samplesheet.csv --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/mag -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

See the [nf-core/mag website documentation](https://nf-co.re/mag/parameters) for more information about pipeline specific parameters.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/mag
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mag releases page](https://github.com/nf-core/mag/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

Additionally, to enable also reproducible results from the individual assembly tools this pipeline provides extra parameters. SPAdes is designed to be deterministic for a given number of threads. To generate reproducible results set the number of cpus with `--spades_fix_cpus` or `--spadeshybrid_fix_cpus`. This will overwrite the number of cpus specified in the `base.config` file and additionally ensure that it is not increased in case of retries for individual samples. MEGAHIT only generates reproducible results when run single-threaded.
You can fix this by using the parameter `--megahit_fix_cpu_1`. In both cases, do not specify the number of cpus for these processes in additional custom config files, this would result in an error.

MetaBAT2 is run by default with a fixed seed within this pipeline, thus producing reproducible results.

Using the BUSCO auto-lineage mode with an internet connection may lead to non-reproducible results, since the databases are frequently updated and automatic lineage selection depends on the version of the database used when running BUSCO.
Therefore, we strongly recommend downloading the required lineage datasets in advance and specifying the lineage to check against.
To ensure reproducibility when using auto-lineage mode, download `all` lineages (see [Databases](#databases)) and provide the download path to `--busco_db`. This will enable offline mode and produce consistent results across runs.

For the taxonomic bin classification with [CAT](https://github.com/dutilh/CAT), when running the pipeline with `--cat_db_generate` the parameter `--save_cat_db` can be used to also save the generated database to allow reproducibility in future runs. Note that when specifying a pre-built database with `--cat_db`, currently the database can not be saved.

The taxonomic classification of bins with GTDB-Tk is not guaranteed to be reproducible, since the placement of bins in the reference tree is non-deterministic. However, the authors of the GTDB-Tk article examined the reproducibility on a set of 100 genomes across 50 trials and did not observe any difference (see [https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)).

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version may be out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

Note, do not change number of CPUs with custom config files for the processes `spades`, `spadeshybrid` or `megahit` when specifying the parameters `--spades_fix_cpus`, `--spadeshybrid_fix_cpus` and `--megahit_fix_cpu_1` respectively.

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

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
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## A note on the ancient DNA subworkflow

nf-core/mag integrates an additional subworkflow to validate ancient DNA _de novo_ assembly:

[Characteristic patterns of ancient DNA (aDNA) damage](<(https://doi.org/10.1073/pnas.0704665104)>), namely DNA fragmentation and cytosine deamination (observed as C-to-T transitions) are typically used to authenticate aDNA sequences. By identifying assembled contigs carrying typical aDNA damages using [PyDamage](https://github.com/maxibor/pydamage), nf-core/mag can report and distinguish ancient contigs from contigs carrying no aDNA damage. Furthermore, to mitigate the effect of aDNA damage on contig sequence assembly, [freebayes](https://github.com/freebayes/freebayes) in combination with [BCFtools](https://github.com/samtools/bcftools) are used to (re)call the variants from the reads aligned to the contigs, and (re)generate contig consensus sequences.

## A note on coverage estimation

In order to run the binning tools included in the pipeline, MAG must first align reads back to the assemblies, and estimate the coverage of each contig.

During the coverage estimation step, these alignments are by default filtered to retain alignments that have a percentage identity of 97% (i.e., of the base pairs that match between the read and the contig, 97% are identical). This value is a good default for short read Illumina data, however for certain long read technologies, the error rates in the reads can be much higher.
For example, older Oxford Nanopore chemistries can have error rates approaching
15% - 20%.

If you are having trouble with the coverage estimation steps (for example, the output depths for each bin are all at or near zero), it may be worth manually adjusting this parameter, if it is appropriate for your data.
You can do this by adjusting the `--longread_percentidentity` and `--shortread_percentidentity` parameters for long reads and short reads, respectively.
For older ONT data, you may wish to look at values of around 85% to improve coverage estimation.

## A note on bin refinement

### Error Reporting

DAS Tool may not always be able to refine bins due to insufficient recovery of enough single-copy genes. In these cases you will get a NOTE such as

```bash
[16/d330a6] NOTE: Process `NFCORE_MAG:MAG:BINNING_REFINEMENT:DASTOOL_DASTOOL (test_minigut_sample2)` terminated with an error exit status (1) -- Error is ignored
```

In this case, DAS Tool has not necessarily failed but was unable to complete the refinement. You will therefore not expect to find any output files in the `GenomeBinning/DASTool/` results directory for that particular sample.

If you are regularly getting such errors, you can try reducing the `--refine_bins_dastool_threshold` value, which will modify the scoring threshold defined in the [DAS Tool publication](https://www.nature.com/articles/s41564-018-0171-1).

## A note on bin filtering

The pipeline offers the ability to filter out bins that fall outside of a certain size in base pairs (`--bin_max_length`, `--bin_min_length`).

This can be useful if you have a set of target organisms that you know approximately the size of the genome for, or if you are looking to filter out small bins that are likely to be contaminants or assembly artifacts.
By removing these bins, you can speed up run time of the pipeline considerably in some cases.

This can also remove 'nonsense' bins of e.g. a single or a collection of very short contigs that can occur with more aggressive binners (e.g. CONCOCT), and can in some cases prevent GUNC [from running correctly](https://github.com/grp-bork/gunc/issues/42#issue-2148763805).

Note that in this context, it is recommended to also set `--min_length_unbinned_contigs` to a suitably high value that corresponds to a reasonable bin size if the `-bin_*_length` parameters are used, so you have useful 'singular' contigs in the unbinned output.

## A note on GTDB having too many files or using too many inodes

The GTDB is very large both in size and by the number of files it contains.
The uncompressed database requires >200k files (or more specifically: inodes), which can be problematic for users with limited storage resources.

One workaround for this is to economize on files/inodes by using a SquashFS image version of the `.tar.gz` GTDB archive, and supply this to the pipeline via a configuration file.

:::warning
This feature is only available with container engines `apptainer` and `singularity`!
:::

To generate your SquashFS image:

1. Install [squashfs-tools](https://github.com/plougher/squashfs-tools), if it is not already on your system
2. Download the GTDB archive either via the full [`.tar.gz` archive](https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/) or the [split database](https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/split_package) version.
3. Convert to a SquashFS image

- Full package (compressed):

  ```bash
  gzip -cd gtdbtk_r220_data.tar.gz | mksquashfs - gtdbtk_r220.squashfs -tar
  ```

- Full package (uncompressed)

  ```bash
  mksquashfs /path/to/database gtdbtk_r220.squashfs
  ```

- Split package (compressed)

  ```bash
  cat gtdbtk_r220_data.tar.gz.part_* | gzip -cd - | mksquashfs - gtdbtk_r220.squashfs -tar
  ```

To use the image in the pipeline:

1. Make an empty directory somewhere on your system
2. Make a custom [Nextflow config](https://nextflow.io/docs/latest/config.html#configuration-file) file with the following settings

   ```nextflow
   process {
       withName: GTDBTK_CLASSIFYWF {
               containerOptions = "-B /<path>/<to>/gtdbtk_r220.squashfs:${params.gtdb_db}:image-src=/"
       }
   }
   ```

3. Supply the empty directory and config to the nextflow run `--gtdb_db /<path>/<to>/<empty_dir>/` to your `nextflow run` command

   ```bash
   nextflow run nf-core/mag -r <version> -profile <profiles> <...> --gtdb_db /<path>/<to>/<empty_dir>/ -c <custom>.config
   ```

:::warning
Make sure to update the paths where indicated, and the GTDB release version if using a more recent one than r220!
:::

:::note
If you have issues with this, you may need to specify a different `image-src=` so it corresponds to the directory structure within your `SquashFS` image.

You can determine this with:

```bash
unsquashfs -l -max-depth 1 -d'' gtdbtk_r220.squashfs
```

And use the resulting output in `image-src=`

```nextflow
process {
    withName: GTDBTK_CLASSIFYWF {
            containerOptions = "-B /<path>/<to>/gtdbtk_r220.squashfs:${params.gtdb_db}:image-src=/<output_from_unsquashfs_ls>"
    }
}
```

Where we update the `image-src` and as above supply the same `/<path>/<to>/<empty_dir>/` path to `--gtdb_db`.
:::

## A note on taxonomic profiling

Generating a taxonomic profile of raw sequencing reads prior to assembly can be highly informative, especially when you are not entirely sure what is in your metagenomic samples.
This can help you identify potential contamination or better understand the taxonomic composition of your samples before proceeding with assembly and binning.

Up until version 4.0.0, this pipeline offered raw read taxonomic profiling using [Centrifuge](https://github.com/centrifugal/centrifuge) and [Kraken2](https://github.com/DerrickWood/kraken2).
This feature was removed in version 5.0.0 to strengthen the pipeline's focus on metagenome assembly and binning.

If you require taxonomic profiling of raw reads, we recommend using [nf-core/taxprofiler](https://nf-co.re/taxprofiler/), which is specifically designed for taxonomic profiling of raw reads and supports a wide range of tools for this purpose.

## BIgMAG compatibility

With the parameter `--generate_bigmag_file` a module will be triggered to generate a file that contains the output from all of the bin-quality tools that can be uploaded to the [BIgMAG](https://github.com/jeffe107/BIgMAG) dashboard for visualising and evaluating MAGs.
Please note that generating this file requires the parameters `--run_busco`, `--run_gunc` and `--run_checkm2`, and GTDBTk should be executed (i.e., not skipped).
The file `bigmag_summary.tsv` located at `GenomeBinning/BIgMAG` is the only file needed to run the BIgMAG dashboard.
