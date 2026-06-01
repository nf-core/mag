# New to nf-core/mag?

**nf-core/mag** is a bioinformatics best-practice analysis pipeline for the assembly, binning, and annotation of metagenomes.

As a domain-agnostic pipeline, it is highly flexible and thus very powerful for metagenomic _de novo_ assembly.
This flexibility allows researchers to reuse the pipeline across many different studies.
However, due to the extensive number of steps and options offered to the user, the pipeline can be very overwhelming for newcomers.

This page is aimed providing context and explicit guidance first-time users for their particular use case, orienting them towards useful or relevant options and settings.

:::warning
Some of the guidance here can be subjective based on our own experiences on our own research projects.
Please always consider carefully whether the guidance applies to your use case, for example reviewing the literature for publications on a similar topic as your own.
:::

This page is split into four specific sections:

- [**Brief introduction to nf-core/mag**](#brief-introduction-to-nf-coremag): this section gives a brief overview of what the entire nf-core/mag pipeline does by default
- [**Defaults and tool selection**](#defaults-and-tool-selection): discusses default steps of the pipeline that are executed in more detail, with additional tool suggestions and parameters to explore
- [**Input types**](#input-types): discusses suitable input data configurations depending on your input data files
- [**Domain and research specific guidance**](domain-and-research-specific-guidance): discusses research-domain specific options (e.g. when targeting specific organisms, or types of DNA)

## Brief introduction to nf-core/mag

The primary aim of nf-core/mag is to generate metagenomic assembled genomes (MAG), perform quality control, and evaluate the quality of each MAG.

A run of the pipeline without any customisation or additional parameters will:

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
  - By binning with `MetaBat2`, `MaxBin2`, `CONCOCT`, `COMEBin`, `MetaBinner`, `SemiBin2`
- Post-binning tasks
  - Quality control bins with `QUAST` and `BUSCO`
  - Annotate bins with `PROKKA` (for bacteria and archaea) and `MetaEuk` (for Eukaryotes)
  - Taxonomic assign bins with `GTDB-Tk` (for bacteria and archaea)

The pipeline will download database files for you where necessary.
All other functionality must be activated using dedicated [parameters](https://nf-co.re/mag/parameters).

The rest of this section provides recommendations and advice for tool selection in these default sections of the pipeline.

## Defaults and tool selection

The pipeline's defaults are designed to perform comprehensive quality control, support all appropriate assemblers and binners, and employ broadly applicable bin evaluation tools.

The defaults prioritize exploring all options to find the optimal tool combination rather than optimizing for speed or resource efficiency.
The following sections provide non-exhaustive ideas for adapting tool choices if desired or needed.

### Preprocessing

Preprocessing of short and long reads does _not_ include host depletion by default. If samples were derived from a host organism (for example, mouse or human), removing host-derived data can improve runtime and assembly quality.

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

The following table summarizes available assemblers:

| Assembler      | Input              | Comment                                                                                                                          |
| -------------- | ------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| `MEGAHIT`      | Short reads        | Fast and memory-efficient; produces competitive assemblies though occasionally with higher misassembly rates compared to SPAdes. |
| `SPAdes`       | Short reads        | Computationally demanding but produces high-quality assemblies with lower misassembly rates; slower than MEGAHIT.                |
| `SPAdesHybrid` | Short & long reads | Slow and memory-intensive; leverages both read types for improved assembly accuracy.                                             |
| `FLYE`         | Long reads         | Slow and memory-intensive; suitable for long-read assemblies but not optimized for speed.                                        |
| `MetaDBG`      | Long reads         | Fast and memory-efficient alternative to FLYE for long-read assembly.                                                            |

When both short and long reads are available, consider running `SPAdesHybrid` and/or **long-read assembly with FLYE or MetaDBG**.

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

The following table summarizes available binners:

| Binner      | Comment                                                                                                                                                                   |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MetaBat2`  | Unsupervised probabilistic binner combining sequence composition and differential coverage; performs well in multi-sample mode and is widely used in ensemble pi`pelines. |
| `MaxBin2`   | Uses Expectation-Maximization with tetranucleotide frequency and single-copy marker genes; particularly effective with multiple sa`mples.                                 |
| `CONCOCT`   | Unsupervised Gaussian mixture model clustering with strong performance in multi-sample binning; frequently complemented by other binners in pi`pelines.                   |
| `COMEBin`   | Deep learning-based binner optimized for contrastive multi-view learning; shows strong performance on hybrid and long-read as`semblies.                                   |
| `Metabinner | Machine learning approach combining contig composition and coverage; designed for improved accuracy in complex me`tagenomes.                                              |
| `Semibin2`  | Semi-supervised binner using pre-trained models and marker genes; performs competitively on diverse metagenome types including challenging datasets.                      |

All binners currently run exclusively with CPUs.
GPU-based execution should accelerate several binners considerably.

> ![NOTE]
> CONCOCT and COMEBin typically have long run times, and not recommended for time sensitive projects, although can have better quality bins in some cases.

**Recommendation**: run as many binners as computationally feasible.

### Should I refine my bins?

Bin refinement tools aims to cross-compare outputs from all binners, and select the highest quality 'version' of the same bin from across the binners.

Bin refinement can improve genome recovery by consolidating outputs from all binners and selecting the "best" result using DAS Tool ([Song and Thomas, 2017](https://doi.org/10.1093/bioinformatics/btx086)).

Bin refinement in nf-core/mag is optional and not enabled by default.

**Recommendation**: select on a per-project basis.

### Which bin quality control tool should I use?

Bin quality controls checks for quality of the assembled bins against a range of criteria such as contamination and completeness, for example as described in the MIMAG reporting standard (
[Bowers et al. 2017](https://doi.org/10.1038/nbt.3893)).

Bin quality control is performed by default in nf-core/mag using BUSCO, which evaluates both prokaryotes and eukaryotes based on marker genes.

Alternatives include **CheckM** (marker-gene-based, prokaryote-only) and **CheckM2** (machine learning-based, prokaryote-only).
CheckM is well established and commonly used, but the newer CheckM2 version has an more recent database and accurate evaluation.

Changing the default is typically driven by the need for comparability with other studies.

Additional chimerism checks with GUNC can be enabled if desired.

**Recommendation**: select on a per-project basis - BUSCO for projects with a broader taxonomic target, CheckM/CheckM2 for prokaryote only targets.

### Which bin taxonomic assignment tool should I use?

Taxonomic assignment involves comparing bins for similarity to known genomes and other MAGs.

**GTDBTk** classifies bins using specific marker genes and yields GTDB-based taxonomies; this approach requires bins of at least medium quality for accuracy. GTDB only supports bacteria and archaea.
**CAT**, by contrast, uses all detectable genes to assign NCBI-based taxonomies. CAT only supports microbes.

The choice between GTDBTk and CAT depends on the desired taxonomy framework or the completeness of bins.

> ![WARNING]
> All bin taxonomic assignment tools in nf-core mag requires very large reference databases (10s-100GBs)!

**Recommendation**: select on a per-project basis.

## Input types

<!-- TODO JAmes continue from here -->

### Short reads assembly

Short read assembly is performed by default if you have specified paths to short-read FASTQ files in the [input samplesheet](../usage.md#input-specifications), and you do not turn off running MEGAHIT and SPAdes using the relevant `--skip_*` [parameters](https://nf-co.re/mag/parameters/).

### Long reads assembly

Long read assembly is performed by default if you have specified long-read FASTQ files in the [input samplesheet](../usage.md#input-specifications), and you do not turn off running Flye and MetaDBG using the relevant `--skip_*` [parameters](https://nf-co.re/mag/parameters/).

### Hybrid assembly (short and long reads)

Hybrid-read assembly is performed by default if you have specified _both_ short- and long-read FASTQ files in the [input samplesheet](../usage.md#input-specifications), and you do not turn off running SPAdesHbyrid using the relevant `--skip_*` [parameters](https://nf-co.re/mag/parameters/).

### Pre-assembled input (binning only)

## Domain and research specific guidance

### Virus Identification

**Virus identification and ancient DNA analysis are not enabled by default.** Virus identification is worth considering if viruses may be relevant to your research question.

### Eukaryotic MAGs

**Annotation of MAGs for eukaryotic genes is not enabled by default.** You must select a database with `--metaeuk_mmseqs_db` to activate this.

### Ancient DNA

**Ancient DNA analysis are not enabled by default.** Ancient DNA processing involves specialized tools and analysis steps to ensure best-practice handling of that particular data type.
