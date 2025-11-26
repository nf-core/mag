# nf-core/mag: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v5.3.0dev [unreleased]

### `Added`

- [#905](https://github.com/nf-core/mag/pull/905) - Add nf-test snapshot for `test_assembly_input` profile (by @dialvarezs)
- [#930](https://github.com/nf-core/mag/pull/930) - Add binner SemiBin2 (by @d4straub)

### `Changed`

- [#932](https://github.com/nf-core/mag/pull/932) - Replaced usages of deprecated `Channel()` with `channel()` and fix other LSP warnings (by @dialvarezs)
- [#937](https://github.com/nf-core/mag/pull/937) - Updated to nf-core 3.5.1 `TEMPLATE` (by @dialvarezs)
- [#938](https://github.com/nf-core/mag/pull/938) - Update nf-core modules (by @dialvarezs)

### `Fixed`

- [#894](https://github.com/nf-core/mag/pull/894) - Fix read order in metaSPAdes to allow co-assembly of paired-end data of multiple samples (reported by @maartenciers, fix by @jfy133 with contributions from @prototaxites, @d4straub and @dialvarezs)
- [#927](https://github.com/nf-core/mag/pull/927) - MetaBinner now succeeds when no contigs are too short or all are binned (reported by @MicroSeq, fix by @d4straub)
- [#929](https://github.com/nf-core/mag/pull/929) - Allow the domain_classification.R script to run with any assembler, not just Megahit or Spades (reported by @MicroSeq, fix by @prototaxites)

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| bcftools | 1.21             | 1.22        |
| fastp    | 0.24.0           | 1.0.1       |
| geNomad  | 1.11.1           | 1.11.2      |
| metamdbg | 1.1              | 1.2         |
| mmseqs   | 17.b804f         | 18.8cc5c    |
| nf-core  |                  | 3.5.1       |
| samtools | 1.21             | 1.22.1      |
| SemiBin2 |                  | 2.2.0       |

### `Deprecated`

## v5.2.0 - Puce Pangolin [2025-11-07]

### `Added`

- [#842](https://github.com/nf-core/mag/pull/842) - Add support for running multiple binQC tools in one run using dedicated `--run_busco`, `--run_checkm`, and `--run_checkm2` parameters (by @harper357, with contributions from @dialvarezs, @prototaxites and @jfy133)
- [#881](https://github.com/nf-core/mag/pull/881) - Add binner MetaBinner (by @d4straub, insprired by @HeshamAlmessady & @AlphaSquad)

### `Changed`

- [#842](https://github.com/nf-core/mag/pull/842) - Change `bin_summary.tsv` format for improved clarity and more comprehensiveness (by @harper357, with contributions from @dialvarezs, @prototaxites and @jfy133)
  - Now will include columns from all bin QC tools executed in a given run (i.e., all/any of BUSCO, CheckM and CheckM2)
  - Adds suffixes to all columns (`_<toolname>`) to distinguish which column comes from which tool

### `Fixed`

- [#896](https://github.com/nf-core/mag/pull/896) - Remove obsolete execution command from README (by @dialvarezs)
- [#907](https://github.com/nf-core/mag/pull/907) - Include refined bins from all binners in the `DASTool/bins` output folder (by @AlexHoratio)
- [#911](https://github.com/nf-core/mag/pull/911) - Ensure column order is consistent when generating depth summaries to prevent swapped results on merged depth summary (by @dialvarezs)
- [#912](https://github.com/nf-core/mag/pull/912) - Fix validation of multiple sequencing platforms when using `binning_map_mode = "all"` (reported by @mjfi2sb3, fix by @dialvarezs)
- [#921](https://github.com/nf-core/mag/pull/921) - Fix publishing of BUSCO files (reported by @joao1980, fix by @dialvarezs)

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| MetaBinner |                  | 1.4.4-0     |

### `Deprecated`

- [#842](https://github.com/nf-core/mag/pull/842) - Remove `--binqc_tool` (by @harper357, with contributions from @dialvarezs, @prototaxites and @jfy133)

## v5.1.0 - Platinum Pudu [2025-10-27]

### `Added`

- [#873](https://github.com/nf-core/mag/pull/873) - Document usage of `longread_percentidentity` and `shortread_percentidentity` and set the value of `longread_percentidentity` in the `test_full` profile to 85 (by @prototaxites)
- [#875](https://github.com/nf-core/mag/pull/875) - Add binner COMEBin (by @d4straub)

### `Changed`

- [#878](https://github.com/nf-core/mag/pull/878) - Refine test_full config with optimised resource usage for AWS release megatests (by @jfy133)
- [#880](https://github.com/nf-core/mag/pull/880) - Updated to nf-core 3.4.1 `TEMPLATE` (by @jfy133)

### `Fixed`

- [#878](https://github.com/nf-core/mag/pull/878) - Fix METASPADES process not receiving the correct number of cpus from the fix CPUs parameter (by @jfy133)
- [#885](https://github.com/nf-core/mag/pull/885) - Fix typo in long-read assembly mode selection (reported by @feixiang1209, fix by @jfy133)
- [#888](https://github.com/nf-core/mag/pull/888) - Only error if all bins are size filtered if bins have actually been generated (reported by @hkaspersen, fix by @prototaxites)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| nf-core | 3.3.2            | 3.4.1       |
| COMEBin |                  | 1.0.4       |

### `Deprecated`

## v5.0.0 - Green Squirrel [2025-09-30]

### `Added`

- [#718](https://github.com/nf-core/mag/pull/718) - Add support for independent long-read metagenomic assembly (requested by @ljmesi and many others, added by @muabnezor)
- [#718](https://github.com/nf-core/mag/pull/718) - Added metaMDBG and (meta)Flye as long read assemblers (added by @muabnezor)
- [#718](https://github.com/nf-core/mag/pull/718) - Added host removal for long reads using minimap2 as aligner (added by @muabnezor)
- [#827](https://github.com/nf-core/mag/pull/827) - Added nf-test CI testing for all test profiles (added by @jfy133)
- [#829](https://github.com/nf-core/mag/pull/829) - Add `--skip_shortread_qc` and `--skip_longread_qc` params for skipping certain default preprocessing steps (added by @erikrikarddaniel)
- [#846](https://github.com/nf-core/mag/pull/846) - Improve documentation of `group` samplesheet column (added by @vinisalazar)
- [#855](https://github.com/nf-core/mag/pull/855) - Add basic nf-tests for test_longreadonly, test_longreadonly_alternatives, test_hybrid and test_assembly_input (added by @dialvarezs)
- [#864](https://github.com/nf-core/mag/pull/864) - Add `--gtdbtk_skip_aniscreen` to disable fast classification of genomes by ANI using skani in GTDB-Tk (by @jfy133 and @prototaxites).

### `Changed`

- [#718](https://github.com/nf-core/mag/pull/718) - Refactored all assembly steps into subworkflows (added by @muabnezor)
- [#799](https://github.com/nf-core/mag/pull/799) - Add `--cat_classify_unbinned`, to enable taxonomic classification of unbinned contigs using CAT (requested by @amizeranschi, added by @dialvarezs)
- [#799](https://github.com/nf-core/mag/pull/799) - Upgraded to latest version of CAT_pack modules (requested by @maxibor, added by @dialvarezs)
- [#811](https://github.com/nf-core/mag/pull/811) - Update util modules, and remove aria2 module to replace with native Nextflow downloading of CheckM database (by @dialvarezs)
- [#816](https://github.com/nf-core/mag/pull/816) - Removed all leftover references to conda 'defaults' channel (by @jfy133)
- [#823](https://github.com/nf-core/mag/pull/823) - Updated to nf-core 3.3.1`TEMPLATE` (by @jfy133 )
- [#827](https://github.com/nf-core/mag/pull/827) - Updated to nf-core 3.3.2`TEMPLATE` (by @dialvarezs)
- [#841](https://github.com/nf-core/mag/pull/841) - MultiQC config updated to support CheckM, CheckM2, and GTDB-Tk (by @harper357)
- [#844](https://github.com/nf-core/mag/pull/843) - Change loading mechanism of internal PhiX/Lambda databases to improve Dev UX when schema building (by @jfy133)
- [#851](https://github.com/nf-core/mag/pull/851) - Improve structure of local modules and subworkflows (by @dialvarezs)
- [#853](https://github.com/nf-core/mag/pull/853) - Update nf-core modules and subworkflows (by @dialvarezs)
- [#856](https://github.com/nf-core/mag/pull/856) - Update more nf-core modules (by @dialvarezs)

### `Fixed`

- [#843](https://github.com/nf-core/mag/pull/843) - Fixed issue with large format Bowtie2 index files not being emitted from index module (reported by Nick Eckersley, fix by @jfy133)
- [#847](https://github.com/nf-core/mag/pull/847) - Allow the BBNorm process to use only 0.8 of the memory allocated to the task to stop if from oversubscribing memory (reported by and fix by @erikrikarddaniel)
- [#850](https://github.com/nf-core/mag/pull/850) - Fixed some modules of the GDTBTk subworkflow not being represented in version lists (fix by @jfy133)
- [#852](https://github.com/nf-core/mag/pull/852) - Fixed version reporting by ensure all modules are represented in final version.yml for MultiQC (by @jfy133)
- [#854](https://github.com/nf-core/mag/pull/854) - Update porechop/abi to a patched version to prevent duplicated read names (reported by @palec87, fix by @jfy133)
- [#858](https://github.com/nf-core/mag/pull/858) - Fix a single parameter validation failure reporting errors for all parameters by updated nf-schema to 2.5.1 (reported by @Pranjal-Bioinfo, fix by @nvnieuwk and @jfy133)
- [#864](https://github.com/nf-core/mag/pull/864) - Fix missing multi-threading of MetaEuk easypredict (reported by @OlivierCoen, fix by @prototaxites).

### `Dependencies`

| Tool         | Previous version | New version |
| ------------ | ---------------- | ----------- |
| bcftools     | 1.17             | 1.21        |
| BUSCO        | 5.8.3            | 6.0.0       |
| CAT          | 5.2.3            | 6.0.1       |
| centrifuge   | 1.0.4.1          | 1.0.4.2     |
| dastool      | 1.1.6            | 1.1.7       |
| nanolyse     | 1.41.6           | 1.44.1      |
| fastp        | 0.23.4           | 0.24.0      |
| flye         |                  | 2.9.5       |
| Freebayes    | 1.3.6            | 1.3.10      |
| geNomad      | 1.5.2            | 1.11.0      |
| GTDB-Tk      | 2.4.0            | 2.5.2       |
| metabat2     | 2.15             | 2.17        |
| metamdbg     |                  | 1.0         |
| minimap2     |                  | 2.29        |
| mmseqs2      | 14.7e284         | 17.b804f    |
| samtools     |                  | 1.21        |
| nf-core      | 3.2.0            | 3.3.2       |
| pydamage     | 0.7.0            | 1.0.0       |
| seqtk        | 1.3              | 1.4         |
| porechop_abi | 0.5.0            | 0.5.0post1  |
| NanoPlot     | 1.44.1           | 1.46.1      |

### `Deprecated`

- [#799](https://github.comf/nf-core/mag/pull/799) - Removed `--cat_official_taxonomy` in favour of `--cat_allow_unofficial_lineages` to control CAT's use of unofficial lineages (added by @dialvarezs)
- [#825](https://github.com/nf-core/mag/pull/825) - Removed `--centrifuge_db`, `--kraken2_db`, `--krona_db` and `--skip_krona` parameters following the removal of taxonomic profiling functionality. See nf-core/taxprofiler for replacement (added by @dialvarezs)
- [#851](https://github.com/nf-core/mag/pull/851) - Remove `POOL_READ_*` local modules in favor of nf-core cat/fastq (by @dialvarezs)
- [#855](https://github.com/nf-core/mag/pull/855) - Remove test_adapterremoval, test_ancient_dna, test_bbnorm, test_busco_auto, test_host_rm, test_hybrid_host_rm, test_binrefinement, test_concoct and test_longread profiles (added by @dialvarezs)
- [#864](https://github.com/nf-core/mag/pull/864) - Remove `--gtdb_mash` due to dropping of support by GTDBTk itself (by @prototaxites and @jfy133)

## v4.0.0 - Blue Huemul [2025-05-22]

### `Added`

- [#730](https://github.com/nf-core/mag/pull/730) - Added `--busco_db_lineage` to allow specifying a specific lineage for BUSCO database (added by @jfy133, @dialvarezs)
- [#730](https://github.com/nf-core/mag/pull/730) - Added a new documentation section on database setup (by @jfy133, @dialvarezs)
- [#784](https://github.com/nf-core/mag/pull/784) - Added `--bin_min_size` and `--bin_max_size` parameters to filter out bins based on size (requested by @maxibor, @alexhbnr, added by @jfy133, @prototaxites)
- [#793](https://github.com/nf-core/mag/pull/793) - Document use of a SquashFS image with `--gtdb_db`, useful for limited inode infrastructure (by @muniheart)
- [#805](https://github.com/nf-core/mag/pull/793) - Add support for fastp's `--trim_poly_g` option (by @jfy133)

### `Changed`

- [#730](https://github.com/nf-core/mag/pull/730) - Migrate from local BUSCO module to nf-core one, updating version (by @dialvarezs)
- [#730](https://github.com/nf-core/mag/pull/730) - Use BUSCO database from nf-core test datasets (by @dialvarezs)
- [#788](https://github.com/nf-core/mag/pull/788) - Tweak method of loading GTDB database in GTDBTK_CLASSIFYWF for more stability (reported by @alexhbnr, fix by @jfy133)
- [#800](https://github.com/nf-core/mag/pull/800) - Default branch is now set to `main` (by @jfy133 and @mirpedrol)
- [#801](https://github.com/nf-core/mag/pull/800) - Increase CheckM memory requests to match [recommended requirements](https://github.com/nf-core/mag/pull/800) (by @jfy133)

### `Fixed`

- [#789](https://github.com/nf-core/mag/pull/789) - Improve `--bowtie2_mode` description to clarify default settings (reported by @IceGreb, fix by @jfy133)
- [#798](https://github.com/nf-core/mag/pull/798) - Fix overly strict database validation for `--metauk_db` and improve documentation (reported by @ruqse, fix by @jfy133)
- [#804](https://github.com/nf-core/mag/pull/804) - Fix broken memory specification for FASTQC (reported by @jmichaelegana, fix by @awgymer & @jfy133)

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| BUSCO    | 5.4.3            | 5.8.3       |
| csvtk    |                  | 0.31.0      |
| nextflow | 24.04.2          | 25.04.2     |

### `Deprecated`

- [#730](https://github.com/nf-core/mag/pull/730) - Remove `--busco_auto_lineage_prok` due to update and simplified usage of BUSCO (added by @jfy133, @dialvarezs)

## v3.4.0 - Green Gecko [2025-04-04]

### `Added`

- [#745](https://github.com/nf-core/mag/pull/745) - Added pipeline parameter `spades_downstreaminput` to use contigs instead of scaffolds (by @Pranjal-Bioinfo, @jfy133, @GallVp & @sateeshperi).
- [#745](https://github.com/nf-core/mag/pull/745) - Added `trimmomatic` as an additional pre-processing tool (by @Pranjal-Bioinfo, @jfy133, @GallVp & @sateeshperi).
- [#745](https://github.com/nf-core/mag/pull/745) - Added parameters for `concoct/cut_up_fasta.py` including `bin_concoct_chunksize`, `bin_concoct_overlap` and `bin_concoct_donotconcatlast` (by @Pranjal-Bioinfo, @jfy133, @GallVp & @sateeshperi).
- [#777](https://github.com/nf-core/mag/pull/777) - Improved input validation through additional JSON keywords and error messages (by @agusinac)
- [#784](https://github.com/nf-core/mag/pull/784) - Add ability to filter out bins based on bin size for downstream using `--bin_min_size`, `--bin_max_size` (reported by @maxibor & @alexhbnr, fix by @jfy133)

### `Changed`

- [#774](https://github.com/nf-core/mag/pull/774) - Update CheckM2 to v1.1.0 and default database (by @dialvarezs).

### `Fixed`

- [#726](https://github.com/nf-core/mag/pull/726) - Fix formatting errors to follow Nextflow best practice (by @dialvarezs).
- [#769](https://github.com/nf-core/mag/pull/769) - Fix megahit not emitting correct filenames due to suboptimal arguments ordering (reported and fix by @IceGreb)
- [#771](https://github.com/nf-core/mag/pull/771) - Fix misspecified checkm2 database parameter check (reported by @dpelegri and fix by @jfy133)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| CheckM2 | 1.0.2            | 1.1.0       |
| SPAdes  | 4.0.0            | 4.1.0       |

### `Deprecated`

## v3.3.1 - [2025-02-13]

### `Added`

- [#758](https://github.com/nf-core/mag/pull/758) - Added new diagram in metro-map style (by @jfy133, @prototaxites, @d4straub)

### `Changed`

- [#731](https://github.com/nf-core/mag/pull/747) - Updated to nf-core 3.1.2 `TEMPLATE` (by @jfy133)
- [#755](https://github.com/nf-core/mag/pull/755) - Updated to nf-core 3.2.0 `TEMPLATE` (by @jfy133)

### `Fixed`

- [#748](https://github.com/nf-core/mag/pull/748) - Fix broken phix reference channel when skipping phix removal (reported by @amizeranschi, fix by @muabnezor)
- [#752](https://github.com/nf-core/mag/pull/752) - Fix QUAST results not being displayed when skipping certain steps (reported by @amizeranschi, fix by @jfy133)
- [#753](https://github.com/nf-core/mag/pull/753) - Fix iGenomes reference support for host removal reference genome (reported by @Thomieh73, fix by @jfy133)
- [#759](https://github.com/nf-core/mag/pull/758) - Fixed parameters that allow both files or directories to not error with directories, and general file input validation improvements (reported by @mjfi2sb3, fix by @jfy133)

### `Deprecated`

## v3.3.0 - Red Reindeer [2024-12-19]

### `Added`

- [#692](https://github.com/nf-core/mag/pull/692) - Added Nanoq as optional longread filtering tool (added by @muabnezor)
- [#692](https://github.com/nf-core/mag/pull/692) - Added chopper as optional longread filtering tool and/or phage lambda removal tool (added by @muabnezor)
- [#707](https://github.com/nf-core/mag/pull/707) - Make Bin QC a subworkflow (added by @dialvarezs)
- [#707](https://github.com/nf-core/mag/pull/707) - Added CheckM2 as an alternative bin completeness and QC tool (added by @dialvarezs)
- [#708](https://github.com/nf-core/mag/pull/708) - Added `--exclude_unbins_from_postbinning` parameter to exclude unbinned contigs from post-binning processes, speeding up Prokka in some cases (added by @dialvarezs)
- [#732](https://github.com/nf-core/mag/pull/732) - Added support for Prokka's compliance mode with `--prokka_with_compliance --prokka_compliance_centre <xyz>` (reported by @audy and @Thomieh73, added by @jfy133)

### `Changed`

- [#731](https://github.com/nf-core/mag/pull/731) - Updated to nf-core 3.1.0 `TEMPLATE` (by @jfy133)

### `Fixed`

- [#707](https://github.com/nf-core/mag/pull/708) - Fixed channel passed as GUNC input (added by @dialvarezs)
- [#724](https://github.com/nf-core/mag/pull/724) - Fix quoting in `utils_nfcore_mag_pipeline/main.nf` (added by @dialvarezs)
- [#716](https://github.com/nf-core/mag/pull/692) - Make short read processing a subworkflow (added by @muabnezor)
- [#708](https://github.com/nf-core/mag/pull/708) - Fixed channel passed as GUNC input (added by @dialvarezs)
- [#729](https://github.com/nf-core/mag/pull/729) - Fixed misspecified multi-FASTQ input for single-end data in MEGAHIT (reported by John Richards, fix by @jfy133)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| CheckM  | 1.2.1            | 1.2.3       |
| CheckM2 |                  | 1.0.2       |
| chopper |                  | 0.9.0       |
| GUNC    | 1.0.5            | 1.0.6       |
| nanoq   |                  | 0.10.0      |

### `Deprecated`

## v3.2.1 [2024-10-30]

### `Added`

### `Changed`

### `Fixed`

- [#707](https://github.com/nf-core/mag/pull/674) - Fix missing space resulting in malformed args for MEGAHIT (reported by @d4straub, fix by @jfy133)

### `Dependencies`

### `Deprecated`

## v3.2.0 - Salmon Salmon [2024-10-27]

### `Added`

- [#674](https://github.com/nf-core/mag/pull/674) - Added `--longread_adaptertrimming_tool` Where user can chose between porechop_abi (default) and porechop (added by @muabnezor)

### `Changed`

- [#674](https://github.com/nf-core/mag/pull/674) - Changed to porechop-abi as default adapter trimming tool for long reads. User can still use porechop if preferred (added by @muabnezor)
- [#666](https://github.com/nf-core/mag/pull/666) - Update SPAdes to version 4.0.0, replace both METASPADES and MEGAHIT with official nf-core modules (requested by @elsherbini, fix by @jfy133)
- [#666](https://github.com/nf-core/mag/pull/666) - Update URLs to GTDB database downloads due to server move (reported by @Jokendo-collab, fix by @jfy133)
- [#695](https://github.com/nf-core/mag/pull/695) - Updated to nf-core 3.0.2 `TEMPLATE` (by @jfy133)
- [#695](https://github.com/nf-core/mag/pull/695) - Switch more stable Zenodo link for CheckM data (by @jfy133)

### `Fixed`

- [#674](https://github.com/nf-core/mag/pull/674) - Make longread preprocessing a subworkflow (added by @muabnezor)
- [#674](https://github.com/nf-core/mag/pull/674) - Add porechop and filtlong logs to multiqc (added by @muabnezor)
- [#674](https://github.com/nf-core/mag/pull/674) - Change local filtlong module to the official nf-core/filtlong module (added by @muabnezor)
- [#690](https://github.com/nf-core/mag/pull/690) - MaxBin2 now using the abundance information from different samples rather than an average (reported by @uel3 and fixed by @d4straub)
- [#698](https://github.com/nf-core/mag/pull/698) - Updated prodigal module to not pick up input symlinks for compression causing pigz errors (reported by @zackhenny, fix by @jfy133 )

### `Dependencies`

| Tool         | Previous version | New version |
| ------------ | ---------------- | ----------- |
| Porechop_ABI |                  | 0.5.0       |
| Filtlong     | 0.2.0            | 0.2.1       |
| SPAdes       | 3.15.3           | 4.0.0       |

### `Deprecated`

## v3.1.0 - Pink Panda [2024-10-04]

### `Added`

- [#665](https://github.com/nf-core/mag/pull/648) - Add support for supplying pre-made bowtie host reference index (requested by @simone-pignotti, added by @jfy133)
- [#670](https://github.com/nf-core/mag/pull/670) - Added `--gtdbtk_pplacer_useram` to run GTDBTk in memory mode rather than write to disk (requested by @harper357, fixed by @jfy133)

### `Changed`

- [#664](https://github.com/nf-core/mag/pull/664) - Update GTDBTk to latest version, with updated column names, update GTDB to release 220 (by @dialvarezs)
- [#676](https://github.com/nf-core/mag/pull/676) - Added exit code 12 to valid SPAdes retry codes, due to OOM errors from spades-hammer (reported by @bawee, fix by @jfy133)

### `Fixed`

- [#667](https://github.com/nf-core/mag/pull/667) - Fix pipeline crashing if only CONCOCT selected during binning (reported and fixed by @jfy133)
- [#670](https://github.com/nf-core/mag/pull/670) - Re-add missing GTDBTk parameters into GTDBTk module (reported by harper357, fixed by @jfy133)
- [#672](https://github.com/nf-core/mag/pull/673) - Fix GTDB-Tk per-sample TSV files not being published in output directory (reported by @jhayer, fix by @jfy133)

### `Dependencies`

| Tool   | Previous version | New version |
| ------ | ---------------- | ----------- |
| GTDBTk | 2.3.2            | 2.4.0       |

### `Deprecated`

- [#670](https://github.com/nf-core/mag/pull/670) - Deprecated `--gtdbtk_pplacer_scratch` due to unintuitive usage (reported by harper357, fixed by @jfy133)

## v3.0.3 [2024-08-27]

### `Added`

### `Changed`

### `Fixed`

- [#648](https://github.com/nf-core/mag/pull/648) - Fix sample ID/assembly ID check failure when no IDs match (reported by @zackhenny, fix by @prototaxites)
- [#646](https://github.com/nf-core/mag/pull/646) - GTDB-Tk directory input now creates a value channel so it runs for all entries to the process and not just the first (reported by @amizeranschi, fix by @prototaxites).
- [#639](https://github.com/nf-core/mag/pull/639) - Fix pipeline failure when a sample produces only a single bin (fix by @d-callan)
- [#651](https://github.com/nf-core/mag/pull/651) - Replace base container for bash only modules to reduce number of containers in pipeline (reported and fixed by @harper357)
- [#652](https://github.com/nf-core/mag/pull/652) - Fix documentation typo in using user-defined assembly parameters (reported and fixed by @amizeranschi)
- [#653](https://github.com/nf-core/mag/pull/653) - Fix overwriting of per-bin 'raw' GUNC RUN output files (multi-bin summary tables not affected) (reported by @zackhenny and fixed by @jfy133)

### `Dependencies`

### `Deprecated`

## v3.0.2 [2024-07-04]

### `Added`

### `Changed`

- [#633](https://github.com/nf-core/mag/pull/633/) - Changed BUSCO to use offline mode when the database is specified by the user (reported by @ChristophKnapp and many others, fix by @jfy133)
- [#632](https://github.com/nf-core/mag/pull/632) - Use default NanoLyse log of just removed reads rather than custom (by @jfy133)

### `Fixed`

- [#630](https://github.com/nf-core/mag/pull/630) - Fix CONCOCT empty bins killing the pipeline, and allow for true multithreading again (removing OPENBLAS loop) (reported by @maxibor, fix by @maxibor and @jfy133)

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| Porechop | 0.2.3_seqan2.1.1 | 0.2.4       |
| NanoPlot | 1.26.3           | 1.41.6      |
| NanoLyse | 1.1.0            | 1.2.0       |

### `Deprecated`

## v3.0.1 [2024-06-10]

### `Added`

### `Changed`

- [#625](https://github.com/nf-core/mag/pull/625) - Updated link to geNomad database for downloading (reported by @amizeranschi, fix by @jfy133)

### `Fixed`

- [#618](https://github.com/nf-core/mag/pull/618) - Fix CENTRIFUGE mkfifo failures by using work directory /tmp (reported by @skrakau, fix by @jfy133)

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| Centrifuge | 1.0.4_beta       | 1.0.4.1     |

### `Deprecated`

## v3.0.0 - Magenta Magpie [2024-05-13]

### `Added`

- [#615](https://github.com/nf-core/mag/pull/615) - Add new logo (by @jfy133)

### `Changed`

- [#599](https://github.com/nf-core/mag/pull/599) - Update to nf-core v2.13.1 `TEMPLATE` (by @jfy133)
- [#614](https://github.com/nf-core/mag/pull/614) - Update to nf-core v2.14.1 `TEMPLATE` (by @jfy133)

### `Fixed`

- [#606](https://github.com/nf-core/mag/pull/606) - Prevent pipeline crash when premade mashdb given to or no alignments found with GTDB-TK_CLASSIFYWF (reported by @cedwardson4, fix by @jfy133)

### `Dependencies`

### `Deprecated`

- [#599](https://github.com/nf-core/mag/pull/599) - Direct reads input (`--input 'sample_{R1,R2}.fastq.gz'`) is no longer supported, all input must come via samplesheets (by @jfy133)

## v2.5.4 [2024-02-12]

### `Added`

### `Changed`

- [#581](https://github.com/nf-core/mag/pull/581) - Added explicit licence text to headers of all custom scripts (reported by @FriederikeHanssen and @maxibor, fix by @jfy133)
- [#602](https://github.com/nf-core/mag/pull/602) - Co-binning when using aDNA mode now enabled (added by @maxibor)

### `Fixed`

- [#583](https://github.com/nf-core/mag/pull/583) - Fix GTDB database input when directory supplied (fix by @jfy133)

### `Dependencies`

### `Deprecated`

## v2.5.3 [2024-02-05]

### `Added`

### `Changed`

- [#575](https://github.com/nf-core/mag/pull/575) - Deactivated MetaSPAdes, Centrifuge, and GTDB in test_full profile due to some container incompatibilities in nf-core megatest AWS configurations (by @jfy133)

### `Fixed`

- [#574](https://github.com/nf-core/mag/pull/574) - Fix wrong channel going to BIN_SUMMARY (fix by @maxibor)

### `Dependencies`

### `Deprecated`

## v2.5.2 [2024-02-02]

### `Added`

- [#562](https://github.com/nf-core/mag/pull/562) - Add CAT summary into the global bin_summary (by @maxibor)
- [#565](https://github.com/nf-core/mag/pull/565) - Add warning of empty GTDB-TK results if no contigs pass completeness filter (by @jfy133 and @maxibor)

### `Changed`

- [#563](https://github.com/nf-core/mag/pull/562) - Update to nf-core v2.12 `TEMPLATE` (by @CarsonJM)
- [#566](https://github.com/nf-core/mag/pull/566) - More logical ordering of MultiQC sections (assembly and bin sections go together respectively) (fix by @jfy133)

### `Fixed`

- [#548](https://github.com/nf-core/mag/pull/548) - Fixes to (reported by @maxibor, @PPpissar, @muniheart, @llborcard, fix by @maxibor)
  - GTDBK-TK execution
  - CAT/QUAST/DEPTH bin summary file name collisions
  - BUSCO database parsing
  - Correct CAT name files
- [#558](https://github.com/nf-core/mag/pull/558) - Fix bug in run merging when dealing with single end data (reported by @roberta-davidson, fix by @jfy133)

### `Dependencies`

### `Deprecated`

## v2.5.1 [2023-11-17]

### `Added`

### `Changed`

### `Fixed`

- [#489](https://github.com/nf-core/mag/pull/489) - Fix file name collision clashes for CHECKM, CAT, GTDBTK, and QUAST (reported by @tillenglert and @maxibor, fix by @maxibor)
- [#533](https://github.com/nf-core/mag/pull/533) - Fix glob pattern for publishing MetaBAT2 bins in results (reported by @patriciatran, fix by @jfy133)
- [#535](https://github.com/nf-core/mag/pull/535) - Fix input validation pattern to again allow direct FASTQ input (reported by @lennijusten, @emnilsson, fix by @jfy133, @d4straub, @mahesh-panchal, @nvnieuwk)

### `Dependencies`

| Tool | Previous version | New version |
| ---- | ---------------- | ----------- |
| CAT  | 4.6              | 5.2.3       |

### `Deprecated`

- [#536](https://github.com/nf-core/mag/pull/536) - Remove custom function with native Nextflow for checking file extension (reported by @d4straub, fix by @jfy133)

## 2.5.0 - Aquamarine Kangaroo - [2023-10-10]

### `Added`

- [#504](https://github.com/nf-core/mag/pull/504) - New parameters `--busco_db`, `--kraken2_db`, and `--centrifuge_db` now support directory input of a pre-uncompressed database archive directory (by @gregorysprenger).

### `Changed`

- [#511](https://github.com/nf-core/mag/pull/511) - Update to nf-core 2.10 `TEMPLATE` (by @jfy133)
- [#504](https://github.com/nf-core/mag/pull/504) - `--save_busco_reference` is now replaced by `--save_busco_db` (by @gregorysprenger).

### `Fixed`

- [#514](https://github.com/nf-core/mag/pull/514) - Fix missing CONCOCT files in downstream output (reported by @maxibor, fix by @jfy133)
- [#515](https://github.com/nf-core/mag/pull/515) - Fix overwriting of GUNC output directories when running with domain classification (reported by @maxibor, fix by @jfy133)
- [#516](https://github.com/nf-core/mag/pull/516) - Fix edge-case bug where MEGAHIT re-uses previous work directory on resume and fails (reported by @husensofteng, fix by @prototaxites)
- [#520](https://github.com/nf-core/mag/pull/520) - Fix missing Tiara output files (fix by @jfy133)
- [#522](https://github.com/nf-core/mag/pull/522) - Fix 'nulls' in depth plot PNG files (fix by @jfy133)

### `Dependencies`

### `Deprecated`

- [#504](https://github.com/nf-core/mag/pull/504) - `--busco_reference`, `--busco_download_path`, `--save_busco_reference` parameters have been deprecated and replaced with new parameters (by @gregorysprenger).

## v2.4.0 - Grey Hammerhead [2023-09-26]

### `Added`

- [#497](https://github.com/nf-core/mag/pull/497) - Adds support for pointing at a local db for krona, using the parameter `--krona_db` (by @willros).
- [#395](https://github.com/nf-core/mag/pull/395) - Adds support for fast domain-level classification of bins using Tiara, to allow bins to be separated into eukaryotic and prokaryotic-specific processes.
- [#422](https://github.com/nf-core/mag/pull/422) - Adds support for normalization of read depth with BBNorm (added by @erikrikarddaniel and @fabianegli)
- [#439](https://github.com/nf-core/mag/pull/439) - Adds ability to enter the pipeline at the binning stage by providing a CSV of pre-computed assemblies (by @prototaxites)
- [#459](https://github.com/nf-core/mag/pull/459) - Adds ability to skip damage correction step in the ancient DNA workflow and just run pyDamage (by @jfy133)
- [#364](https://github.com/nf-core/mag/pull/364) - Adds geNomad nf-core modules for identifying viruses in assemblies (by @PhilPalmer and @CarsonJM)
- [#481](https://github.com/nf-core/mag/pull/481) - Adds MetaEuk for annotation of eukaryotic MAGs, and MMSeqs2 to enable downloading databases for MetaEuk (by @prototaxites)
- [#437](https://github.com/nf-core/mag/pull/429) - `--gtdb_db` also now supports directory input of an pre-uncompressed GTDB archive directory (reported by @alneberg, fix by @jfy133)
- [#494](https://github.com/nf-core/mag/pull/494) - Adds support for saving the BAM files from Bowtie2 mapping of input reads back to assembly (fix by @jfy133)

### `Changed`

- [#428](https://github.com/nf-core/mag/pull/428) [#467](https://github.com/nf-core/mag/pull/467) - Update to nf-core 2.8, 2.9 `TEMPLATE` (by @jfy133)
- [#429](https://github.com/nf-core/mag/pull/429) - Replaced hardcoded CheckM database auto-download URL to a parameter (reported by @erikrikarddaniel, fix by @jfy133)
- [#441](https://github.com/nf-core/mag/pull/441) - Deactivated CONCOCT in AWS 'full test' due to very long runtime (fix by @jfy133).
- [#442](https://github.com/nf-core/mag/pull/442) - Remove warning when BUSCO finds no genes in bins, as this can be expected in some datasets (reported by @Lumimar, fix by @jfy133).
- [#444](https://github.com/nf-core/mag/pull/444) - Moved BUSCO bash code to script (by @jfy133)
- [#477](https://github.com/nf-core/mag/pull/477) - `--gtdb` parameter is split into `--skip_gtdbtk` and `--gtdb_db` to allow finer control over GTDB database retrieval (fix by @jfy133)
- [#500](https://github.com/nf-core/mag/pull/500) - Temporarily disabled downstream processing of both refined and raw bins due to bug (by @jfy133)

### `Fixed`

- [#496](https://github.com/nf-core/mag/pull/496) - Fix help text for paramters `--bowtie2_mode`, `spades_options` and `megahit_options` (by @willros)
- [#400](https://github.com/nf-core/mag/pull/400) - Fix duplicated Zenodo badge in README (by @jfy133)
- [#406](https://github.com/nf-core/mag/pull/406) - Fix CheckM database always downloading, regardless if CheckM is selected (by @jfy133)
- [#419](https://github.com/nf-core/mag/pull/419) - Fix bug with busco_clean parameter, where it is always activated (by @prototaxites)
- [#426](https://github.com/nf-core/mag/pull/426) - Fixed typo in help text for parameters `--host_genome` and `--host_fasta` (by @tillenglert)
- [#434](https://github.com/nf-core/mag/pull/434) - Fix location of samplesheet for AWS full tests (reported by @Lfulcrum, fix by @jfy133)
- [#438](https://github.com/nf-core/mag/pull/438) - Fixed version inconsistency between conda and containers for GTDBTK_CLASSIFYWF (by @jfy133)
- [#439](https://github.com/nf-core/mag/pull/445) - Fix bug in assembly input (by @prototaxites)
- [#447](https://github.com/nf-core/mag/pull/447) - Remove `default: None` from parameter schema (by @drpatelh)
- [#449](https://github.com/nf-core/mag/pull/447) - Fix results file overwriting in Ancient DNA workflow (reported by @alexhbnr, fix by @jfy133)
- [#470](https://github.com/nf-core/mag/pull/470) - Fix binning preparation from running even when binning was requested to be skipped (reported by @prototaxites, fix by @jfy133)
- [#480](https://github.com/nf-core/mag/pull/480) - Improved `-resume` reliability through better meta map preservation (reported by @prototaxites, fix by @jfy133)
- [#493](https://github.com/nf-core/mag/pull/493) - Update `METABAT2` nf-core module so that it reduced the number of unnecessary file moves, enabling virtual filesystems (fix by @adamrtalbot)
- [#500](https://github.com/nf-core/mag/pull/500) - Fix MaxBin2 bins not being saved in results directly properly (reported by @Perugolate, fix by @jfy133)

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| BCFtools | 1.16             | 1.17        |
| SAMtools | 1.16.1           | 1.17        |
| fastp    | 0.23.2           | 0.23.4      |
| MultiQC  | 1.14             | 1.15        |

## v2.3.2 - [2023-06-23]

### `Fixed`

- [#461](https://github.com/nf-core/mag/pull/461) - Fix full-size AWS test profile paths (by @jfy133)
- [#461](https://github.com/nf-core/mag/pull/461) - Fix pyDamage results being overwritten (reported by @alexhbnr, fix by @jfy133)

## v2.3.1 - [2023-06-19]

### `Fixed`

- [#458](https://github.com/nf-core/mag/pull/458) - Correct the major issue in ancient DNA workflow of binning refinement being performed on uncorrected contigs instead of aDNA consensus recalled contigs (issue [#449](https://github.com/nf-core/mag/issues/449))
- [#451](https://github.com/nf-core/mag/pull/451) - Fix results file overwriting in Ancient DNA workflow (reported by @alexhbnr, fix by @jfy133, and integrated by @maxibor in [#458](https://github.com/nf-core/mag/pull/458) )

## v2.3.0 - Red Cow [2023-03-02]

### `Added`

- [#350](https://github.com/nf-core/mag/pull/350) - Adds support for CheckM as alternative bin completeness and QC tool (added by @jfy133 and @skrakau)
- [#353](https://github.com/nf-core/mag/pull/353) - Added the busco_clean parameter to optionally clean each BUSCO directory after a successful (by @prototaxites)
- [#361](https://github.com/nf-core/mag/pull/361) - Added the skip_clipping parameter to skip read preprocessing with fastp or adapterremoval. Running the pipeline with skip_clipping, keep_phix and without specifying a host genome or fasta file skips the FASTQC_TRIMMED process (by @prototaxites)
- [#365](https://github.com/nf-core/mag/pull/365) - Added CONCOCT as an additional (optional) binning tool (by @jfy133)
- [#366](https://github.com/nf-core/mag/pull/366) - Added CAT_SUMMARISE process and cat_official_taxonomy parameter (by @prototaxites)
- [#372](https://github.com/nf-core/mag/pull/372) - Allow CAT_DB to take an extracted database as well as a tar.gz file (by @prototaxites).
- [#380](https://github.com/nf-core/mag/pull/380) - Added support for saving processed reads (clipped, host removed etc.) to results directory (by @jfy133)
- [#394](https://github.com/nf-core/mag/pull/394) - Added GUNC for additional chimeric bin/contamination QC (added by @jfy133)

### `Changed`

- [#340](https://github.com/nf-core/mag/pull/340),[#368](https://github.com/nf-core/mag/pull/368),[#373](https://github.com/nf-core/mag/pull/373) - Update to nf-core 2.7.2 `TEMPLATE` (by @jfy133, @d4straub, @skrakau)
- [#373](https://github.com/nf-core/mag/pull/373) - Removed parameter `--enable_conda`. Updated local modules to new conda syntax and updated nf-core modules (by @skrakau)
- [#385](https://github.com/nf-core/mag/pull/385) - CAT also now runs on unbinned contigs as well as binned contigs (added by @jfy133)
- [#399](https://github.com/nf-core/mag/pull/399/files) - Removed undocumented BUSCO_PLOT process (previously generated `*.busco_figure.png` plots unsuitable for metagenomics) (by @skrakau).
- [#416](https://github.com/nf-core/mag/pull/416) - Use GTDBTK_CLASSIFYWF nf-core module instead of local module (added by @alxndrdiaz)

### `Fixed`

- [#345](https://github.com/nf-core/mag/pull/345) - Bowtie2 mode changed to global alignment for ancient DNA mode (`--very-sensitive` mode) to prevent soft clipping at the end of reads when running in local mode. (by @maxibor)
- [#349](https://github.com/nf-core/mag/pull/349) - Add a warning that pipeline will reset minimum contig size to 1500 specifically MetaBAT2 process, if a user supplies below this threshold. (by @jfy133)
- [#352](https://github.com/nf-core/mag/pull/352) - Escape the case in the BUSCO module that BUSCO can just detect a root lineage but is not able to find any marker genes (by @alexhbnr)
- [#355](https://github.com/nf-core/mag/pull/355) - Include error code 21 for retrying with higher memory for SPAdes and hybridSPAdes (by @mglubber)

### `Dependencies`

| Tool      | Previous version | New version |
| --------- | ---------------- | ----------- |
| BUSCO     | 5.1.0            | 5.4.3       |
| BCFtools  | 1.14             | 1.16        |
| Freebayes | 1.3.5            | 1.3.6       |
| SAMtools  | 1.15             | 1.16.1      |

## v2.2.1 [2022-08-25]

### `Added`

### `Changed`

### `Fixed`

- [#328](https://github.com/nf-core/mag/pull/328) - Fix too many symbolic links issue in local convert_depths module (reported by @ChristophKnapp and fixed by @apeltzer, @jfy133)
- [#329](https://github.com/nf-core/mag/pull/329) - Each sample now gets it's own result directory for PyDamage analysis and filter (reported and fixed by @maxibor)

### `Dependencies`

## v2.2.0 - Golden Mammoth [2022-06-14]

### `Added`

- [#263](https://github.com/nf-core/mag/pull/263) - Restructure binning subworkflow in preparation for aDNA workflow and extended binning
- [#247](https://github.com/nf-core/mag/pull/247) - Add ancient DNA subworkflow
- [#263](https://github.com/nf-core/mag/pull/263) - Add MaxBin2 as second contig binning tool
- [#285](https://github.com/nf-core/mag/pull/285) - Add AdapterRemoval2 as an alternative read trimmer
- [#291](https://github.com/nf-core/mag/pull/291) - Add DAS Tool for bin refinement
- [#319](https://github.com/nf-core/mag/pull/319) - Activate pipeline-specific institutional nf-core/configs

### `Changed`

- [#269](https://github.com/nf-core/mag/pull/269),[#283](https://github.com/nf-core/mag/pull/283),[#289](https://github.com/nf-core/mag/pull/289),[#302](https://github.com/nf-core/mag/pull/302) - Update to nf-core 2.4 `TEMPLATE`
- [#286](https://github.com/nf-core/mag/pull/286) - Cite our publication instead of the preprint
- [#291](https://github.com/nf-core/mag/pull/291), [#299](https://github.com/nf-core/mag/pull/299) - Add extra results folder `GenomeBinning/depths/contigs` for `[assembler]-[sample/group]-depth.txt.gz`, and `GenomeBinning/depths/bins` for `bin_depths_summary.tsv` and `[assembler]-[binner]-[sample/group]-binDepths.heatmap.png`
- [#315](https://github.com/nf-core/mag/pull/315) - Replace base container for standard shell tools to fix problems with running on Google Cloud

### `Fixed`

- [#290](https://github.com/nf-core/mag/pull/290) - Fix caching of binning input
- [#305](https://github.com/nf-core/mag/pull/305) - Add missing Bowtie2 version for process `BOWTIE2_PHIX_REMOVAL_ALIGN` to `software_versions.yml`
- [#307](https://github.com/nf-core/mag/pull/307) - Fix retrieval of GTDB-Tk version (note about newer version caused error in `CUSTOM_DUMPSOFTWAREVERSIONS`)
- [#309](https://github.com/nf-core/mag/pull/309) - Fix publishing of BUSCO `busco_downloads/` folder, i.e. publish only when `--save_busco_reference` is specified
- [#321](https://github.com/nf-core/mag/pull/321) - Fix parameter processing in `BOWTIE2_REMOVAL_ALIGN` (which was erroneously for `BOWTIE2_PHIX_REMOVAL_ALIGN`)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| fastp   | 0.20.1           | 0.23.2      |
| MultiQC | 1.11             | 1.12        |

## v2.1.1 [2021-11-25]

### `Added`

- [#240](https://github.com/nf-core/mag/pull/240) - Add prodigal to predict protein-coding genes for assemblies.
- [#241](https://github.com/nf-core/mag/pull/241) - Add parameter `--skip_prodigal`.
- [#244](https://github.com/nf-core/mag/pull/244) - Add pipeline preprint information.
- [#245](https://github.com/nf-core/mag/pull/245) - Add Prokka to annotate binned genomes.

### `Changed`

- [#249](https://github.com/nf-core/mag/pull/249) - Update workflow overview figure.
- [#258](https://github.com/nf-core/mag/pull/258) - Updated MultiQC 1.9 to 1.11.
- [#260](https://github.com/nf-core/mag/pull/260) - Updated SPAdes 3.13.1 -> 3.15.3, MEGAHIT 1.2.7 -> 1.2.7

### `Fixed`

- [#256](https://github.com/nf-core/mag/pull/256) - Fix `--skip_busco`.
- [#236](https://github.com/nf-core/mag/pull/236) - Fix large assemblies (> 4 billion nucleotides in length).
- [#254](https://github.com/nf-core/mag/pull/254) - Fix MetaBAT2 error with nextflow version 21.10.x (21.04.03 is the latest functional version for nf-core/mag 2.1.0).
- [#255](https://github.com/nf-core/mag/pull/255) - Update gtdbtk conda channel.
- [#258](https://github.com/nf-core/mag/pull/258) - FastP results are now in MultiQC.

## v2.1.0 - Black Zebra [2021-07-29]

### `Added`

- [#212](https://github.com/nf-core/mag/pull/212), [#214](https://github.com/nf-core/mag/pull/214) - Add bin abundance estimation based on median sequencing depths of corresponding contigs (results are written to `results/GenomeBinning/bin_depths_summary.tsv` and `results/GenomeBinning/bin_summary.tsv`) [#197](https://github.com/nf-core/mag/issues/197).
- [#214](https://github.com/nf-core/mag/pull/214) - Add generation of (clustered) heat maps with bin abundances across samples (using centered log-ratios)
- [#217](https://github.com/nf-core/mag/pull/217) - Publish genes predicted with Prodigal within BUSCO run (written to `results/GenomeBinning/QC/BUSCO/[assembler]-[bin]_prodigal.gff`).

### `Changed`

- [#218](https://github.com/nf-core/mag/pull/218) - Update to nf-core 2.0.1 `TEMPLATE` (DSL2)

### `Fixed`

- [#226](https://github.com/nf-core/mag/pull/226) - Fix handling of `BUSCO` output when run in auto lineage selection mode and selected specific lineage is the same as the generic one.

## v2.0.0 - Silver Swan [2021-06-01]

### `Added`

- [#179](https://github.com/nf-core/mag/pull/179) - Add BUSCO automated lineage selection functionality (new default). The pameter `--busco_auto_lineage_prok` can be used to only consider prokaryotes and the parameter `--busco_download_path` to run BUSCO in `offline` mode.
- [#178](https://github.com/nf-core/mag/pull/178) - Add taxonomic bin classification with `GTDB-Tk` `v1.5.0` (for bins filtered based on `BUSCO` QC metrics).
- [#196](https://github.com/nf-core/mag/pull/196) - Add process for CAT database creation as an alternative to using pre-built databases.

### `Changed`

- [#162](https://github.com/nf-core/mag/pull/162) - Switch to DSL2
- [#162](https://github.com/nf-core/mag/pull/162) - Changed `--input` file format from `TSV` to `CSV` format, requires header now
- [#162](https://github.com/nf-core/mag/pull/162) - Update `README.md`, `docs/usage.md` and `docs/output.md`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `FastP` from version `0.20.0` to `0.20.1`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `Bowtie2` from version `2.3.5` to `2.4.2`
- [#162](https://github.com/nf-core/mag/pull/162) - Update `FastQC` from version `0.11.8` to `0.11.9`
- [#172](https://github.com/nf-core/mag/pull/172) - Compressed discarded MetaBAT2 output files
- [#176](https://github.com/nf-core/mag/pull/176) - Update CAT DB link
- [#179](https://github.com/nf-core/mag/pull/179) - Update `BUSCO` from version `4.1.4` to `5.1.0`
- [#179](https://github.com/nf-core/mag/pull/179) - By default BUSCO now performs automated lineage selection instead of using the bacteria_odb10 lineage as reference. Specific lineage datasets can still be provided via `--busco_reference`.
- [#178](https://github.com/nf-core/mag/pull/178) - Change output file: `results/GenomeBinning/QC/quast_and_busco_summary.tsv` -> `results/GenomeBinning/bin_summary.tsv`, contains GTDB-Tk results as well.
- [#191](https://github.com/nf-core/mag/pull/191) - Update to nf-core 1.14 `TEMPLATE`
- [#193](https://github.com/nf-core/mag/pull/193) - Compress CAT output files [#180](https://github.com/nf-core/mag/issues/180)
- [#198](https://github.com/nf-core/mag/pull/198) - Requires nextflow version `>= 21.04.0`
- [#200](https://github.com/nf-core/mag/pull/200) - Small changes in GitHub Actions tests
- [#203](https://github.com/nf-core/mag/pull/203) - Renamed `fastp` params and improved description in documentation: `--mean_quality` -> `--fastp_qualified_quality`, `--trimming_quality` -> `--fastp_cut_mean_quality`

### `Fixed`

- [#175](https://github.com/nf-core/mag/pull/175) - Fix bug in retrieving the `--max_unbinned_contigs` longest unbinned sequences that are longer than `--min_length_unbinned_contigs` (`split_fasta.py`)
- [#175](https://github.com/nf-core/mag/pull/175) - Improved runtime of `split_fasta.py` in `METABAT2` process (important for large assemblies, e.g. when computing co-assemblies)
- [#194](https://github.com/nf-core/mag/pull/194) - Allow different folder structures for Kraken2 databases containing `*.k2d` files [#187](https://github.com/nf-core/mag/issues/187)
- [#195](https://github.com/nf-core/mag/pull/195) - Fix documentation regarding required compression of input FastQ files [#160](https://github.com/nf-core/mag/issues/160)
- [#196](https://github.com/nf-core/mag/pull/196) - Add process for CAT database creation as solution for problem caused by incompatible `DIAMOND` version used for pre-built `CAT database` and `CAT classification` [#90](https://github.com/nf-core/mag/issues/90), [#188](https://github.com/nf-core/mag/issues/188)

## v1.2.0 - Yellow Squirrel [2021-02-10]

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

### `Removed`

- [#143](https://github.com/nf-core/mag/pull/143) - Change parameter: `--manifest` -> `--input`

## v1.1.2 - Blue Panda [2020-11-24]

### `Changed`

- [#135](https://github.com/nf-core/mag/pull/135) - Update to nf-core 1.12 `TEMPLATE`

### `Fixed`

- [#133](https://github.com/nf-core/mag/pull/133) - Fixed processing of `--input` parameter [#131](https://github.com/nf-core/mag/issues/131)

## v1.1.1 - Lime Owl [2020-11-10]

### `Added`

- [#121](https://github.com/nf-core/mag/pull/121) - Add full-size test
- [#124](https://github.com/nf-core/mag/pull/124) - Add worfklow overview figure to `README`

### `Changed`

- [#123](https://github.com/nf-core/mag/pull/123) - Update to new nf-core 1.11 `TEMPLATE`

### `Fixed`

- [#118](https://github.com/nf-core/mag/pull/118) - Fix `seaborn` to `v0.10.1` to avoid `nanoplot` error
- [#120](https://github.com/nf-core/mag/pull/120) - Fix link to CAT database in help message
- [#124](https://github.com/nf-core/mag/pull/124) - Fix description of `CAT` process in `output.md`

## v1.1.0 - White Elephant [2020-10-06]

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

## v1.0.0 - Purple Corgi [2019-12-20]

Initial release of nf-core/mag, created with the [nf-core](http://nf-co.re/) template.

### `Added`

- short and long reads QC (fastp, porechop, filtlong, fastqc)
- Lambda and PhiX detection and filtering (bowtie2, nanolyse)
- Taxonomic classification of reads (centrifuge, kraken2)
- Short read and hybrid assembly (megahit, metaspades)
- metagenome binning (metabat2)
- QC of bins (busco, quast)
- annotation (cat/bat)
