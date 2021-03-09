# nf-core/mag: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0dev - [date]

### `Added`

### `Changed`

- [#162](https://github.com/nf-core/mag/pull/162) - Switch to DSL2
- [#162](https://github.com/nf-core/mag/pull/162) - Changed `--input` file format from `TSV` to `CSV` format, requires header now
- [#162](https://github.com/nf-core/mag/pull/162) - Update `README.md`, `docs/usage.md` and `docs/output.md`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `FastP` from version `0.20.0` to `0.20.1`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `Bowtie2` from version `2.3.5` to `2.4.2`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `FastQC` from version `0.11.8` to `0.11.9`

## v1.2.0 - 2020/02/10

### `Added`

- [#146](https://github.com/nf-core/mag/pull/146) - Add `--coassemble_group` parameter to allow group-wise co-assembly
- [#146](https://github.com/nf-core/mag/pull/146) - Add `--binning_map_mode` parameter allowing different mapping strategies to compute co-abundances used for binning (`all`, `group`, `own`)
- [#149](https://github.com/nf-core/mag/pull/149) - Add two new parameters to allow custom SPAdes and MEGAHIT options (`--spades_options` and `--megahit_options`)

### `Changed`

- [#141](https://github.com/nf-core/mag/pull/141) - Update to nf-core 1.12.1 `TEMPLATE`
- [#143](https://github.com/nf-core/mag/pull/143) - Manifest file has to be handed over via `--input` parameter now
- [#143](https://github.com/nf-core/mag/pull/143) - Changed format of manifest input file: requires a '.tsv' suffix and additionally contains group ID
- [#143](https://github.com/nf-core/mag/pull/143) - TSV `--input` file allows now also entries containing only short reads
- [#145](https://github.com/nf-core/mag/pull/145) - When using TSV input files, uses sample IDs now for `FastQC` instead of basenames of original read files. Allows non-unique file basenames.

### `Deprecated`

- [#143](https://github.com/nf-core/mag/pull/143) - Change depreciated parameters: `--manifest` -> `--input`

## v1.1.2 - 2020/11/24

### `Changed`

- [#135](https://github.com/nf-core/mag/pull/135) - Update to nf-core 1.12 `TEMPLATE`

### `Fixed`

- [#133](https://github.com/nf-core/mag/pull/133) - Fixed processing of `--input` parameter [#131](https://github.com/nf-core/mag/issues/131)

## v1.1.1 - 2020/11/10

### `Added`

- [#121](https://github.com/nf-core/mag/pull/121) - Add full-size test
- [#124](https://github.com/nf-core/mag/pull/124) - Add worfklow overview figure to `README`

### `Changed`

- [#123](https://github.com/nf-core/mag/pull/123) - Update to new nf-core 1.11 `TEMPLATE`

### `Fixed`

- [#118](https://github.com/nf-core/mag/pull/118) - Fix `seaborn` to `v0.10.1` to avoid `nanoplot` error
- [#120](https://github.com/nf-core/mag/pull/120) - Fix link to CAT database in help message
- [#124](https://github.com/nf-core/mag/pull/124) - Fix description of `CAT` process in `output.md`

## v1.1.0 - 2020/10/06

### `Added`

- [#35](https://github.com/nf-core/mag/pull/35) - Add social preview image
- [#49](https://github.com/nf-core/mag/pull/49) - Add host read removal with `Bowtie 2` and according custom section to `MultiQC`
- [#49](https://github.com/nf-core/mag/pull/49) - Add separate `MultiQC` section for `FastQC` after preprocessing
- [#65](https://github.com/nf-core/mag/pull/65) - Add `MetaBAT2` RNG seed parameter `--metabat_rng_seed` and set the default to 1 which ensures reproducible binning results
- [#65](https://github.com/nf-core/mag/pull/65) - Add parameters `--megahit_fix_cpu_1`, `--spades_fix_cpus` and `--spadeshybrid_fix_cpus` to ensure reproducible results from assembly tools
- [#66](https://github.com/nf-core/mag/pull/66) - Export `depth.txt.gz` into result folder
- [#67](https://github.com/nf-core/mag/pull/67) - Compress assembly files
- [#82](https://github.com/nf-core/mag/pull/82) - Add `nextflow_schema.json`
- [#104](https://github.com/nf-core/mag/pull/104) - Add parameter `--save_busco_reference`

### `Changed`

- [#56](https://github.com/nf-core/mag/pull/56) - Update `MetaBAT2` from `v2.13` to `v2.15`
- [#46](https://github.com/nf-core/mag/pull/46) - Update `MultiQC` from `v1.7` to `v1.9`
- [#88](https://github.com/nf-core/mag/pull/88) - Update to new nf-core 1.10.2 `TEMPLATE`
- [#88](https://github.com/nf-core/mag/pull/88) - `--reads` is now removed, use `--input` instead
- [#101](https://github.com/nf-core/mag/pull/101) - Prevented PhiX alignments from being stored in work directory [#97](https://github.com/nf-core/mag/issues/97)
- [#104](https://github.com/nf-core/mag/pull/104), [#111](https://github.com/nf-core/mag/pull/111) - Update `BUSCO` from `v3.0.2` to `v4.1.4`

### `Fixed`

- [#29](https://github.com/nf-core/mag/pull/29) - Fix `MetaBAT2` binning discards unbinned contigs [#27](https://github.com/nf-core/mag/issues/27)
- [#31](https://github.com/nf-core/mag/pull/31), [#36](https://github.com/nf-core/mag/pull/36), [#76](https://github.com/nf-core/mag/pull/76), [#107](https://github.com/nf-core/mag/pull/107) - Fix links in README
- [#47](https://github.com/nf-core/mag/pull/47) - Fix missing `MultiQC` when `--skip_quast` or `--skip_busco` was specified
- [#49](https://github.com/nf-core/mag/pull/49), [#89](https://github.com/nf-core/mag/pull/89) - Added missing parameters to summary
- [#50](https://github.com/nf-core/mag/pull/50) - Fix missing channels when `--keep_phix` is specified
- [#54](https://github.com/nf-core/mag/pull/54) - Updated links to `minikraken db`
- [#54](https://github.com/nf-core/mag/pull/54) - Fixed `Kraken2` dp preparation: allow different names for compressed archive file and contained folder as for some minikraken dbs
- [#55](https://github.com/nf-core/mag/pull/55) - Fixed channel joining for multiple samples causing `MetaBAT2` error [#32](https://github.com/nf-core/mag/issues/32)
- [#57](https://github.com/nf-core/mag/pull/57) - Fix number of threads used by `MetaBAT2` program `jgi_summarize_bam_contig_depths`
- [#70](https://github.com/nf-core/mag/pull/70) - Fix `SPAdes` memory conversion issue [#61](https://github.com/nf-core/mag/issues/61)
- [#71](https://github.com/nf-core/mag/pull/71) - No more ignoring errors in `SPAdes` assembly
- [#72](https://github.com/nf-core/mag/pull/72) - No more ignoring of `BUSCO` errors
- [#73](https://github.com/nf-core/mag/pull/73), [#75](https://github.com/nf-core/mag/pull/75) - Improved output documentation
- [#96](https://github.com/nf-core/mag/pull/96) - Fix missing bin names in `MultiQC` BUSCO section [#78](https://github.com/nf-core/mag/issues/78)
- [#104](https://github.com/nf-core/mag/pull/104) - Fix `BUSCO` errors causing missing summary output [#77](https://github.com/nf-core/mag/issues/77)

### `Deprecated`

- [#29](https://github.com/nf-core/mag/pull/29) - Change depreciated parameters: `--singleEnd` -> `--single_end`, `--igenomesIgnore` -> `--igenomes_ignore`

## v1.0.0 - 2019/12/20

Initial release of nf-core/mag, created with the [nf-core](http://nf-co.re/) template.

### `Added`

- short and long reads QC (fastp, porechop, filtlong, fastqc)
- Lambda and PhiX detection and filtering (bowtie2, nanolyse)
- Taxonomic classification of reads (centrifuge, kraken2)
- Short read and hybrid assembly (megahit, metaspades)
- metagenome binning (metabat2)
- QC of bins (busco, quast)
- annotation (cat/bat)
