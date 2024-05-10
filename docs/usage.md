# nf-core/mag: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mag/usage](https://nf-co.re/mag/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Input specifications

The input data can be passed to nf-core/mag in two possible ways, either using the `--input` parameter of raw-reads alone or `--input` additionally with `--assembly_input` that specifies pre-built assemblies.

### Samplesheet input file

You can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata. Furthermore when a `run` column is present, the pipeline will also run perform run- or lane-wise concatenation, for cases where you may have a sample or library sequenced with the same sequencing configuration across multiple runs. The optional run merging happens after short read QC (adapter clipping, host/PhiX removal etc.), and prior to normalisation, taxonomic profiling, and assembly.

At a minimum CSV file should contain the following columns:

`sample,group,short_reads_1,short_reads_2,long_reads`

The path to `long_reads` and `short_reads_2` is optional. Valid examples could look like the following:

```csv title="samplesheet.csv"
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
```

or

```csv title="samplesheet.csv"
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1.fastq.gz,,
sample2,0,data/sample2.fastq.gz,,
```

or to additionally to perform run merging of two runs of sample1:

```csv title="samplesheet.csv"
sample,run,group,short_reads_1,short_reads_2,long_reads
sample1,1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample1,2,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample2,0,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz
sample3,1,0,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
```

Please note the following requirements:

- a minimum 5 of comma-separated columns
- Valid file extension: `.csv`
- Must contain the header `sample,group,short_reads_1,short_reads_2,long_reads` (where `run` can be optionally added)
- Run IDs must be unique within a multi-run sample. A sample with multiple runs will be automatically concatenated.
- FastQ files must be compressed (`.fastq.gz`, `.fq.gz`)
- `long_reads` can only be provided in combination with paired-end short read data
- Within one samplesheet either only single-end or only paired-end reads can be specified
- If single-end reads are specified, the command line parameter `--single_end` must be specified as well

Again, by default, the group information is only used to compute co-abundances for the binning step, but not for group-wise co-assembly (see the parameter docs for [`--coassemble_group`](https://nf-co.re/mag/parameters#coassemble_group) and [`--binning_map_mode`](https://nf-co.re/mag/parameters#binning_map_mode) for more information about how this group information can be used).

### Supplying pre-computed assemblies

It is also possible to run nf-core/mag on pre-computed assemblies, by supplying a CSV file to the parameter `--assembly_input` in addition to the raw reads supplied to `--input`. Supplying assembly input skips all read pre-processing and assembly, jumping straight to the binning stage of the pipeline.

The assembly CSV file should contain the following columns:

`id,group,assembler,fasta`

Where `id` is the ID of the assembly, group is the assembly/binning group (see samplesheet information section for more details), `assembler` is the assembler used to produce the assembly (one of `MEGAHIT`, `SPAdes`, or `SPAdesHybrid`), and `fasta` is the path to the assembly fasta file. Input fasta files can be compressed or uncompressed, but compressed assemblies will be automatically uncompressed for use within the pipeline. The exact information required for each supplied assembly depends on whether the assemblies provided are single assemblies or group-wise co-assemblies. For the following example `--input` CSV:

```csv title="samplesheet.csv"
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
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

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/mag -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
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

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mag releases page](https://github.com/nf-core/mag/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

Additionally, to enable also reproducible results from the individual assembly tools this pipeline provides extra parameters. SPAdes is designed to be deterministic for a given number of threads. To generate reproducible results set the number of cpus with `--spades_fix_cpus` or `--spadeshybrid_fix_cpus`. This will overwrite the number of cpus specified in the `base.config` file and additionally ensure that it is not increased in case of retries for individual samples. MEGAHIT only generates reproducible results when run single-threaded.
You can fix this by using the prameter `--megahit_fix_cpu_1`. In both cases, do not specify the number of cpus for these processes in additional custom config files, this would result in an error.

MetaBAT2 is run by default with a fixed seed within this pipeline, thus producing reproducible results.

To allow also reproducible bin QC with BUSCO, run BUSCO providing already downloaded lineage datasets (BUSCO will be run using automated lineage selection in offline mode) or provide a specific lineage dataset via `--busco_db` and use the parameter `--save_busco_db`. This may be useful since BUSCO datasets are frequently updated and old versions do not always remain (easily) accessible.

For the taxonomic bin classification with [CAT](https://github.com/dutilh/CAT), when running the pipeline with `--cat_db_generate` the parameter `--save_cat_db` can be used to also save the generated database to allow reproducibility in future runs. Note that when specifying a pre-built database with `--cat_db`, currently the database can not be saved.

When it comes to visualizing taxonomic data using [Krona](https://github.com/marbl/Krona), you have the option to provide a taxonomy file, such as `taxonomy.tab`, using the `--krona_db` parameter. If you don't supply a taxonomy file, Krona is designed to automatically download the required taxonomy data for visualization.

The taxonomic classification of bins with GTDB-Tk is not guaranteed to be reproducible, since the placement of bins in the reference tree is non-deterministic. However, the authors of the GTDB-Tk article examined the reproducibility on a set of 100 genomes across 50 trials and did not observe any difference (see [https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)).

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

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
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
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

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

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

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## A note on the ancient DNA subworkflow

nf-core/mag integrates an additional subworkflow to validate ancient DNA _de novo_ assembly:

[Characteristic patterns of ancient DNA (aDNA) damage](<(https://doi.org/10.1073/pnas.0704665104)>), namely DNA fragmentation and cytosine deamination (observed as C-to-T transitions) are typically used to authenticate aDNA sequences. By identifying assembled contigs carrying typical aDNA damages using [PyDamage](https://github.com/maxibor/pydamage), nf-core/mag can report and distinguish ancient contigs from contigs carrying no aDNA damage. Furthermore, to mitigate the effect of aDNA damage on contig sequence assembly, [freebayes](https://github.com/freebayes/freebayes) in combination with [BCFtools](https://github.com/samtools/bcftools) are used to (re)call the variants from the reads aligned to the contigs, and (re)generate contig consensus sequences.

## A note on bin refinement

### Error Reporting

DAS Tool may not always be able to refine bins due to insufficient recovery of enough single-copy genes. In these cases you will get a NOTE such as

```bash
[16/d330a6] NOTE: Process `NFCORE_MAG:MAG:BINNING_REFINEMENT:DASTOOL_DASTOOL (test_minigut_sample2)` terminated with an error exit status (1) -- Error is ignored
```

In this case, DAS Tool has not necessarily failed but was unable to complete the refinement. You will therefore not expect to find any output files in the `GenomeBinning/DASTool/` results directory for that particular sample.

If you are regularly getting such errors, you can try reducing the `--refine_bins_dastool_threshold` value, which will modify the scoring threshold defined in the [DAS Tool publication](https://www.nature.com/articles/s41564-018-0171-1).
