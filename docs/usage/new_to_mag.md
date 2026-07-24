# New to nf-core/mag?

**nf-core/mag** is a bioinformatics best-practice analysis pipeline for the assembly, binning, and annotation of metagenomes.

As a domain-agnostic pipeline, it is highly flexible and thus very powerful for metagenomic _de novo_ assembly.
This flexibility allows researchers to reuse the pipeline across many different studies.
However, the pipeline can be very overwhelming for newcomers due to the extensive number of steps and options offered.

This page is aimed at providing context and explicit guidance first-time users for their particular use case.
It orients them towards useful or relevant options and settings.

> [!WARNING]
> Some of the guidance here is subjective based on our own experiences on our own research projects.
> Please always consider carefully whether the guidance applies to your use case, for example reviewing the literature for publications on a similar topic as your own.

This page is split into four specific sections:

- [**Brief introduction to nf-core/mag**](#brief-introduction-to-nf-coremag): this section gives a brief overview of what the entire nf-core/mag pipeline does by default.
- [**Defaults and tool selection**](#defaults-and-tool-selection): discusses default steps of the pipeline that are executed in more detail, with additional tool suggestions and parameters to explore.
- [**Input types**](#input-types): discusses suitable input data configurations depending on your input data files.
- [**Domain and research specific guidance**](#domain-and-research-specific-guidance): discusses research-domain specific options (for example, when targeting specific organisms, or types of DNA).

> [!IMPORTANT]
> We highly recommend reading the full sections for each question, as many of the suggestions come with caveats.
> However, we list here a very brief summary of the advice for commonly asked questions:
>
> - **Short read only data?** SPAdes and MEGAHIT will run by default
> - **Long read only data?** Flye and MetaMDBG will run by default
> - **Short and long reads?** All assemblers run by default, including SPAdes in Hybrid mode
> - **Already have assemblies?** Use the `--assembly_input` flag to start from binning
> - **Assemble or co-assemble?** It's complicated! Check our advice below.
> - **Which assembler to use?** Use as many as you can.
> - **Can I polish my assemblies?** No, but it is being considered.
> - **Which binner to use?** Use as many as you can.
> - **Should I refine my bins?** It depends on your context.
> - **Which bin QC tool to use?** CheckM/CheckM2 for prokaryote only projects, BUSCO for projects with wider taxonomic targets.
> - **Which taxonomic assignment tool?** It depends on your context.
> - **I want to find viruses** Explore running with `--run_virus_identification`
> - **I want to identify eukaryotes** Explore running with `--bin_domain_classification` and `--metaeuk_mmseqs_db`
> - **I have ancient DNA** Use `--ancient_dna`

## Brief introduction to nf-core/mag

### What does nf-core/mag do

The primary aim of nf-core/mag is to generate metagenomic assembled genomes (MAGs), perform quality control, and evaluate the quality of each MAG.

A run of nf-core/mag without any customisation or additional parameters will:

- Preprocess input reads to
  - Remove adapters (short reads: `fastp`, long reads: `porechop_ABI`)
  - FASTQ-level quality filtering (long reads: `Filtlong`)
  - Remove Phi X (short-read) or lambda (long-read) sequences
  - Quality control reads (short reads: FastQC, long reads: `Nanoplot`)
- Assemble reads into contigs with
  - `MEGAHIT` and `SPAdes` if short reads are provided
  - `SPAdesHybrid` if short and long reads are provided
  - `Flye` and `metaMDBG` if long reads are provided
- Post-assembly tasks
  - Quality control by `QUAST` and `ALE`
  - Annotate assemblies with `prodigal`
- Group contigs
  - By binning with `MetaBat2`, `MaxBin2`, `CONCOCT`, `COMEBin`, `MetaBinner` and `SemiBin2`
- Post-binning tasks
  - Quality control bins with `QUAST` and `BUSCO`
  - Annotate bins with `PROKKA` (for bacteria and archaea) and `MetaEuk` (for eukaryotes)
  - Taxonomically assign bins with `GTDB-Tk` (for bacteria and archaea)

> [!NOTE]
> If you are new to metagenomic _de novo_ assembly, and while reading this page you do not feel comfortable with the terminology used in this page, we suggest reading some introductory literature before running your analysis.
> A few suggestions are as follows:
>
> - Quince, C., et al. (2017). Shotgun metagenomics, from sampling to analysis. Nature Biotechnology, 35(9), 833–844. [doi:10.1038/nbt.3935](https://doi.org/10.1038/nbt.3935)
> - Goussarov, G., et al. (2022). Introduction to the principles and methods underlying the recovery of metagenome-assembled genomes from metagenomic data. MicrobiologyOpen, 11(3), e1298. [doi:10.1002/mbo3.1298](https://doi.org/10.1002/mbo3.1298)
> - Liu, S., et al. (2025). Analysis of metagenomic data. Nature Reviews. Methods Primers, 5(1), 5. [doi:10.1038/s43586-024-00376-6](https://doi.org/10.1038/s43586-024-00376-6)
> - New, F. N., & Brito, I. L. (2020). What Is Metagenomics Teaching Us, and What Is Missed? Annual Review of Microbiology, 74, 117–135. [doi:10.1146/annurev-micro-012520-072314](https://doi.org/10.1146/annurev-micro-012520-072314)

The pipeline will download database files for you where necessary.
All other functionality must be activated using dedicated [parameters](https://nf-co.re/mag/parameters).

### How do I run nf-core/mag?

Before you start, if you are unfamiliar with executing nf-core pipelines, we highly recommend reading the central [Getting started with nf-core](https://nf-co.re/docs/get_started/nf-core) documentation.

To run nf-core/mag, at a minimum, you will need:

- Nextflow
- A software dependency and compute environment manager (typically Conda, Docker, or Singularity/Apptainer)
- Metagenomic FASTQ files (Illumina or Nanopore sequence data is directly supported)
- A UNIX-based (Linux, Mac OSX etc.) computing infrastructure (for example, laptop, server, high performance computing cluster (HPC), or cloud)
  - Internet connection is needed for automatic download of databases, however you can download these files yourself. See [parameters](https://nf-co.re/mag/parameters/) for more information.
  - You will need a large amount of hard drive space, as some of the the databases (>100 GB per database) and output can be very large

You will then need to prepare a [samplesheet](../usage#samplesheet-input-file) specifying the names and paths to the FASTQ Files.

Finally to execute, you will run a typical Nextflow command, either with [parameters](https://nf-co.re/mag/parameters/) customising the run specified on the command line,

```bash
nextflow run nf-core/mag -r 5.4.2  --input samplesheet.csv --outdir results/ -profile docker --clip_tool adapterremoval --reads_minlength 25 <...>
```

or with parameters specified in a [Nextflow parameter file](https://docs.seqera.io/nextflow/cli#pipeline-parameters).

```bash
nextflow run nf-core/mag -r 5.4.2 -params-file params.json -profile docker
```

Where params.json contains:

```json
{
  "input": "samplesheet.csv",
  "output": "results/",
  "clip_tool": "adapterremoval",
  "reads_minlength": 25
}
```

### What is the output from nf-core/mag?

nf-core/mag produces many different output directories and files, depending on what parameters you set.

Typically, the primary output files you will want from a typical run:

- `multiqc/multiqc_report.html`: For overall run summary statistics
- `GenomeBinning/bin_summary.tsv`: Summary statistics with quality metrics for the bins
- `GenomeBinning/<binner name>/`: The bin FASTA files themselves

You can explore example output from real data on the dedicated [Results](https://nf-co.re/mag/results/) page.

However, you will need to explore output on a per-project basis, as the exact files you will need depend on your question.
You can see detailed descriptions of all output files on the [Output](https://nf-co.re/mag/latest/docs/output/) page.

## Defaults and tool selection

The pipeline's default tool executions are designed to perform comprehensive quality control, support all appropriate assemblers and binners, and employ broadly applicable bin evaluation tools.

The default tool selection prioritises exploring all options to find the optimal tool combination rather than optimising for speed or resource efficiency.
The following sections provide non-exhaustive ideas for adapting tool choices when you desire or need.

### Preprocessing

Preprocessing of short and long reads does _not_ include host depletion by default. If samples were derived from a host organism (for example, mouse, or human), removing host-derived data can improve runtime and assembly quality.

**Recommendation**: Omit or skip preprocessing steps only on a case-by-case basis with documented justification.

### To assemble or co-assemble?

By default, the pipeline assembles each sample separately.

However, _co-assembly_ (pooling data from multiple samples) can be beneficial in some cases by increasing overall sequencing depth and improving comparability.

Conversely, co-assembly can be detrimental when combining unrelated samples, as it increases complexity, drastically raises computational requirements, and promotes chimeric sequences (contigs of mixed origin) ([Hofmeyr et al. 2020](https://doi.org/10.1038/s41598-020-67416-5), [Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)).

**Recommendation**: select on a per-project basis.

### Which assembler to use?

By default all assemblers are run on the provided data.

No assembler performs equally across all datasets, and no tool consistently outperforms others ([Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)).
Therefore, despite the larger computational overhead, testing multiple options can improve results.

However, reducing computational burden and accelerating analysis may require tool selection.

The following table summarises available assemblers:

| Assembler      | Input              | Comment                                                                                                                                                                                 |
| -------------- | ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MEGAHIT`      | Short reads        | Faster and more memory-efficient. Produces competitive assemblies though occasionally with higher misassembly rates compared to SPAdes.                                                 |
| `SPAdes`       | Short reads        | Computationally demanding but produces high-quality assemblies with lower misassembly rates. Slower than MEGAHIT.                                                                       |
| `SPAdesHybrid` | Short & long reads | Slower and more memory-intensive. Leverages both read types for improved assembly accuracy.                                                                                             |
| `Flye`         | Long reads         | Slower and more memory-intensive. Suitable for long-read assemblies but not optimised for speed. In some cases may perform better when you have issues with noisy reads or low coverage |
| `MetaDBG`      | Long reads         | Faster and more memory-efficient alternative to Flye for long-read assembly.                                                                                                            |

When both short and long reads are available, consider running `SPAdesHybrid` and/or **long-read assembly with Flye or MetaDBG**.

With high-depth long reads, long-read assembly typically yields more coherent results ([Agustinho et al. 2024](https://doi.org/10.1038/s41592-024-02262-1)).
Short-read-first assembly performs better with high-depth short reads or low-quality long reads and produces more fragmented but higher-accuracy assemblies ([Overholt et al. 2020](https://doi.org/10.1111/1462-2920.15186), [Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)).

**Recommendation**: run as many assemblers as computationally feasible.

### Can I polish my assemblies?

Polishing of long-read assemblies involves using high-quality short reads to 'repair' mistakes in lower-quality long reads.

Polishing of assemblies with short or long reads is _not_ currently implemented in the pipeline.

For metagenomes, polishing can harm assembly quality by erroneously modifying low-abundance genomes using high-abundance data.
High-quality Nanopore data (10.4) may not benefit substantially from long-read polishing (for example, with `Medaka`).
Polishing long-read assemblies with short-read data might be beneficial but remains debated and is not currently available in the pipeline (for example, `Polypolish` and `Pypolca`).

### What binning tool should I select?

All binners are run by default in nf-core/mag.

All binners use abundance information across one or multiple samples to extract (fragmented) genomes.
This abundance information can be calculated across all samples, a specific sample group, or only the single dataset from which the assembly originates.
Providing more abundance information across samples generally improves binning performance, particularly in multi-sample modes.
However, single-sample binning remains viable when only one dataset is available.

Each binning tool implements different algorithms and approaches, and no binner consistently performs best across all scenarios ([Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)).
Exploring all options can improve results.

The following table summarises available binners:

| Binner       | Comment                                                                                                                                                                  |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `MetaBat2`   | Unsupervised probabilistic binner combining sequence composition and differential coverage. Performs well in multi-sample mode and is widely used in ensemble pipelines. |
| `MaxBin2`    | Uses Expectation-Maximisation with tetranucleotide frequency and single-copy marker genes. Particularly effective with multiple samples.                                 |
| `CONCOCT`    | Unsupervised Gaussian mixture model clustering with strong performance in multi-sample binning. Frequently complemented by other binners in pipelines.                   |
| `COMEBin`    | Deep learning-based binner optimised for contrastive multi-view learning. Shows strong performance on hybrid and long-read assemblies.                                   |
| `Metabinner` | Machine learning approach combining contig composition and coverage. Designed for improved accuracy in complex metagenomes.                                              |
| `Semibin2`   | Semi-supervised binner using pretrained models and marker genes. Performs competitively on diverse metagenome types including challenging datasets.                      |

All binners currently run exclusively with CPUs.
GPU-based execution should accelerate several binners considerably.

> [!NOTE]
> CONCOCT and COMEBin typically have long run times, and are not recommended for time sensitive projects.
> They can produce better quality bins in some cases.

**Recommendation**: run as many binners as computationally feasible.

### Should I refine my bins?

Bin refinement tools aim to cross-compare outputs, and select the highest quality 'version' of the same bin from across all binners.

Bin refinement can improve genome recovery by consolidating outputs from all binners, and selecting the 'best' result using DAS Tool ([Song and Thomas, 2017](https://doi.org/10.1093/bioinformatics/btx086)).

Bin refinement in nf-core/mag is optional and not enabled by default.

**Recommendation**: select on a per-project basis.

### Which bin quality control tool should I use?

Bin quality control checks for quality of the assembled bins against a range of criteria such as contamination and completeness, for example, against the MIMAG reporting standard ([Bowers et al. 2017](https://doi.org/10.1038/nbt.3893)).

Bin quality control is performed by default in nf-core/mag using BUSCO.
BUSCO evaluates both prokaryotes and eukaryotes based on marker genes.

Alternatives include **CheckM** (marker-gene-based, prokaryote-only) and **CheckM2** (machine learning-based, prokaryote-only).
CheckM is well established and commonly used, but the newer CheckM2 version has a more recent database and accurate evaluation.

Changing the default is typically driven by the need for comparability with other studies.

Additional chimerism checks with GUNC can be enabled if desired.

**Recommendation**: select on a per-project basis - BUSCO for projects with a broader taxonomic target, CheckM/CheckM2 for prokaryote only targets.

### Which bin taxonomic assignment tool should I use?

Taxonomic assignment involves comparing bins for similarity to known genomes and other MAGs.

**GTDB-Tk** classifies bins using specific marker genes and yields GTDB-based taxonomies. This approach requires bins of at least medium quality for accuracy. GTDB only supports bacteria and archaea.
**CAT**, by contrast, uses all detectable genes to assign NCBI-based taxonomies. CAT only supports microbes.

The choice between GTDB-Tk and CAT depends on the desired taxonomy framework or the completeness of bins.

> [!WARNING]
> All taxonomic assignment tools for bins in nf-core/mag requires very large reference databases (10-100 GBs)!

**Recommendation**: select on a per-project basis.

### What parameters should I change?

Many of the tools in nf-core/mag can be adjusted using pipeline level parameters (for example, `--min_contig_size` or `--gtdbtk_min_completeness`).

In the vast majority of cases, the defaults to these parameters are set to the same default of the tools themselves.

While these are reasonable defaults, these may not be suitable for all use cases.
These are selected as nf-core/mag is aimed at being a generalist pipeline that can be used in many different research contexts.
You should read the documentation of the specific tools and literature of similar projects to get an idea of which parameters may be useful to adjust for your own project.

**Recommendation**: adjust parameters on a per-project basis.

## Input types

### I only have short reads

Short read assembly is performed by default with MEGAHIT and SPAdes if you include paths to short-read FASTQ files in the [input samplesheet](../usage.md#input-specifications).

You do not need to turn on any assembler.
If you want to turn off running a particular assembler, you can turn it off using the relevant `--skip_*` [parameter](https://nf-co.re/mag/parameters/).

### I only have long reads

Long read assembly is performed by default with Flye and MetaDBG if you have specified long-read FASTQ files in the [input samplesheet](../usage#input-specifications).

You do not need to turn on an assembler.
If you want to turn off running a particular assembler, you can turn it off using the relevant `--skip_*` [parameter](https://nf-co.re/mag/parameters/).

### I have both short and long reads

Hybrid-read assembly is performed by default with SPAdes in hybrid mode if you have specified _both_ short- and long-read FASTQ files in the [input samplesheet](../usage.md#input-specifications), in addition to standalone short-read and long-read assembly.

If you do not want hybrid assembly, you can turn it off using the relevant `--skip_*` [parameter](https://nf-co.re/mag/parameters/) for SPAdes hybrid.

### I already have assemblies, can I just run binning and downstream steps?

nf-core/mag supports 'assembly input', where you can specify pre-assembled contigs in the samplesheet.
This can be useful when you want to use an assembler not currently supported by nf-core/mag, or you wish to reanalyse publicly available assemblies.

For this you must still supply an input samplesheet as usual, however you will skip the assembly step if you provide a second samplesheet to `--assembly_input`.

For more information on how to prepare both samplesheets, see the usage documentation's [Supplying precomputed assemblies](../usage/#supplying-pre-computed-assemblies) section.

## Domain and research specific guidance

### I want to identify viruses in my metagenome

nf-core/mag will not identify viral sequences by default.

If you wish to identify viruses in your metagenome, you need to specify `--run_virus_identification`.

This will execute [`geNomad`](https://github.com/apcamargo/genomad/) to identify contigs that can be assigned to viral genomes.
Note that `geNomad` does not identify viral contigs in bins or MAGs - only raw assemblies.

**Recommendation**: Viral contigs are not screened for by default.
Activate with `--run_virus_identification` to identify viral contigs in assemblies.

### I want to identify eukaryotic MAGs

nf-core/mag will not identify eukaryotic MAGs by default.

If you wish to identify eukaryotic bins your metagenome, you need to run `--bin_domain_classification` and/or specify a `MetaEuk` compatible database.

You can distinguish between eukaryotic and prokaryotic bins if you supply `--bin_domain_classification`, which will execute [Tiara](https://github.com/ibe-uw/tiara).
This will also only send prokaryotic bins to downstream steps that expect only prokaryotic bins.

If a valid prebuilt MetaEuk database name is given to `--metaeuk_mmseqs_db`, or a path supplied to a locally downloaded `MMseqs2` formatted database to `--metaeuk_db`, nf-core/mag will execute [MetaEuk](https://github.com/soedinglab/metaeuk) to predict and annotate eukaryotic genes on bins,

**Recommendation**: Eukaryotic bin identification is not executed by default.
Specify `--bin_domain_classification` to taxonomically classify bins at domain level, or supply a MetaEuk database name or supply a database path to the relevant path to annotate eukaryotic genes in bins.

### I want to assemble Ancient DNA

nf-core/mag will not appropriately analyse ancient DNA sequences by default.

If you have ancient DNA sequences, nf-core/mag has a dedicated sub-workflow to carry out damage pattern authentication and correction.
This does not run by default, and must be activated with `--ancient_dna`.

The pipeline will run `pyDamage` to generate _per-contig_ statistical probabilities of typical ancient DNA deamination miscoding lesion patterns.
If binning is not skipped, you will also receive by default, _per-bin_ averaged damage statistics in the final bin summary file.

If the ancient DNA analysis mode is activated, the pipeline will also by default perform damage correction to remove accidentally incorporated damaged positions in assemblies.

**Recommendation**: Ancient DNA data is not appropriately processed by default. Activate ancient DNA mode with `--ancient_dna` to get damage authentication statistics and correct misincorporated damage in assemblies.
