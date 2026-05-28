# New to nf-core/mag?

**nf-core/mag** is a bioinformatics best-practice analysis pipeline for the assembly, binning, and annotation of metagenomes.

As a domain-agnostic pipeline, it is highly flexible and thus very powerful.
This allows researchers to reuse the pipeline across many different studies.
However, due to the extensive number of steps and options, the pipeline can be very overwhelming for newcomers.

This page is aimed providing context and explicit guidance first-time users to help orient them towards useful options and settings relevant for their particular use case.

:::warning
Some of the guidance here can be subjective based on our own experiences on our own research projects.
Please always consider carefully whether the guidance applies to your usecase.
:::

This page is split into four specific sections:

- [**Brief introduction to nf-core/mag**](#brief-introduction-to-nf-coremag): this section gives a brief overview of what the entire nf-core/mag pipeline does by default
- [**Defaults and tool selection**](#defaults-and-tool-selection): discusses default steps of the pipeline that are executed in more detail, with additional tool suggestions and parameters to explore
- [**Input types**](#input-types): discusses suitable input data configurations depending on your input data files
- [**Domain and research specific guidance**](domain-and-research-specific-guidance): discusses research-domain specific options (e.g. when targeting specific organisms, or types of DNA)

## Brief introduction to nf-core/mag

The primary aim of nf-core/mag is to generate metagenomic assembled genomes (MAG), perform quality control, and evaluate the quality of each MAG.

**The pipeline's defaults are designed to perform comprehensive quality control, support all appropriate assemblers and binners, and employ broadly applicable bin evaluation tools.** The defaults prioritize exploring all options to find the optimal tool combination rather than optimizing for speed or resource efficiency. The following sections provide non-exhaustive ideas for adapting tool choices if desired or needed.

## Defaults and tool selection

### Preprocessing

**Omit or skip preprocessing steps only on a case-by-case basis with documented justification.**

Preprocessing of short and long reads does **not** include host depletion by default. If samples were derived from a host organism (for example, mouse or human), removing host-derived data can improve runtime and assembly quality.

### Assembly

**By default, the pipeline assembles each sample separately.** However, **co-assembly** (pooling data from multiple samples) can be beneficial by increasing overall sequencing depth and improving comparability. Conversely, co-assembly can be detrimental when combining unrelated samples, as it increases complexity, drastically raises computational requirements, and promotes chimeric sequences (contigs of mixed origin) ([Hofmeyr et al. 2020](https://doi.org/10.1038/s41598-020-67416-5), [Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)).

**All assemblers are run on the provided data.** No assembler performs equally across all datasets, and no tool consistently outperforms others ([Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)). Therefore, testing multiple options can improve results. However, reducing computational burden and accelerating analysis may require tool selection. The following table summarizes available assemblers:

| Assembler     | Input              | Comment                                                                                                                          |
| ------------- | ------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| MEGAHIT       | Short reads        | Fast and memory-efficient; produces competitive assemblies though occasionally with higher misassembly rates compared to SPAdes. |
| SPAdes        | Short reads        | Computationally demanding but produces high-quality assemblies with lower misassembly rates; slower than MEGAHIT.                |
| SPAdes Hybrid | Short & long reads | Slow and memory-intensive; leverages both read types for improved assembly accuracy.                                             |
| FLYE          | Long reads         | Slow and memory-intensive; suitable for long-read assemblies but not optimized for speed.                                        |
| MetaDBG       | Long reads         | Fast and memory-efficient alternative to FLYE for long-read assembly.                                                            |

When both short and long reads are available, consider running **SPAdes Hybrid** and/or **long-read assembly with FLYE or MetaDBG**. With high-depth long reads, long-read assembly typically yields more coherent results ([Agustinho et al. 2024](https://doi.org/10.1038/s41592-024-02262-1)). Short-read-first assembly performs better with high-depth short reads or low-quality long reads and produces more fragmented but higher-accuracy assemblies ([Overholt et al. 2020](https://doi.org/10.1111/1462-2920.15186), [Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)).

**No polishing of assemblies with short or long reads is implemented in the pipeline.** For metagenomes, polishing can harm assembly quality by erroneously modifying low-abundance genomes using high-abundance data. High-quality Nanopore data (10.4) may not benefit substantially from long-read polishing (for example, with Medaka). Polishing long-read assemblies with short-read data might be beneficial but remains debated and is not currently available in the pipeline (for example, Polypolish and Pypolca).

### Binning

**All binners use abundance information across one or multiple samples to extract (fragmented) genomes.** This abundance information can be calculated across all samples, a specific sample group, or only the single dataset from which the assembly originates. Providing more abundance information across samples generally improves binning performance, particularly in multi-sample modes. However, single-sample binning remains viable when only one dataset is available.

**All binners are run by default.** These tools implement different algorithms and approaches, and no binner consistently performs best across all scenarios ([Meyer et al. 2022](https://doi.org/10.1038/s41592-022-01431-4)). Exploring all options can improve results. The following table summarizes available binners:

| Binner     | Comment                                                                                                                                                                  |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| MetaBat2   | Unsupervised probabilistic binner combining sequence composition and differential coverage; performs well in multi-sample mode and is widely used in ensemble pipelines. |
| MaxBin2    | Uses Expectation-Maximization with tetranucleotide frequency and single-copy marker genes; particularly effective with multiple samples.                                 |
| CONCOCT    | Unsupervised Gaussian mixture model clustering with strong performance in multi-sample binning; frequently complemented by other binners in pipelines.                   |
| COMEBin    | Deep learning-based binner optimized for contrastive multi-view learning; shows strong performance on hybrid and long-read assemblies.                                   |
| Metabinner | Machine learning approach combining contig composition and coverage; designed for improved accuracy in complex metagenomes.                                              |
| Semibin2   | Semi-supervised binner using pre-trained models and marker genes; performs competitively on diverse metagenome types including challenging datasets.                     |

All binners currently run exclusively with CPUs. GPU-based execution should accelerate several binners considerably.

### Bin Refinement (opt-in)

**Bin refinement can improve genome recovery by consolidating outputs from all binners and selecting the "best" result using DAS Tool** ([Song and Thomas, 2017](https://doi.org/10.1093/bioinformatics/btx086)). Bin refinement is optional and not enabled by default.

### Bin Quality Control and Classification

**Bin quality control is performed by default using BUSCO, which evaluates both prokaryotes and eukaryotes based on marker genes.** Alternatives include **CheckM** (marker-gene-based, prokaryote-only) and **CheckM2** (machine learning-based, prokaryote-only).Changing the default is typically driven by the need for comparability with other studies.

Chimerism checks with GUNC can be enabled if desired.

**Taxonomic classification requires large reference databases.** When supplied, **GTDBTk** classifies bins using specific marker genes and yields GTDB-based taxonomies; this approach requires bins of at least medium quality for accuracy. **CAT**, by contrast, uses all detectable genes to assign NCBI-based taxonomies. The choice between GTDBTk and CAT depends on the desired taxonomy framework or the completeness of bins.

## Input types

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
