# Resource usage guidance

**nf-core/mag** is a bioinformatics best-practice analysis pipeline for the assembly, binning, and annotation of metagenomes.

Due to the nature of metagenomic samples, nf-core/mag can in some cases require large amounts of computational resources to execute.

In this page we provide information and guidance on how nf-core/mag works by default, and in some cases how to optimise computational resource usage.

## Multi-sample or bin executing modules

By default, nf-core/mag aims to make the runtime as efficiently as possible on HPC and similar infrastructure through generating more but shorter running jobs.

In most cases, all pre-assembly steps typically run one job per input FASTQ file.
After assembly, most steps run on a per assembly or contig.
After binning, each step runs on a per-bin basis, and so forth.

There are a few exceptions to this.
In the following cases, depending on the options you set in the pipeline, you may see some steps running as a single job with multiple samples.

- Multiple samples processed in one job:
  - If `--coassemble_group`: all reads across all samples will be pooled into one assembly job (e.g. MEGAHIT, metaSPAdes)
  - If `--gtdbtk_single_job`: all bins across all samples will be pooled into a single GTDBTk classification job

## Default resource requests

The following table lists the default resources that Nextflow will request from a machine on it's first execution attempt for a given module.
It summaries information recorded in the file `conf/base.conf` and the module files themselves `modules`.

To customise these resource requests, see the central [nf-core instructions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources).

Please also note the following:

- Not all modules in this table will run in every pipeline run, nor that the tool will necessarily actually _use_ the requested amount.
- On certain module execution failures (such as out of memory), Nextflow will try and resubmit with increased resources with the equation `<resource value> * task.attempt`, up to a certain number of retries. See `conf/base.conf` for more information.
- Some modules have a dedicate flag to fix the number of CPUs for reproducibility reasons (e.g. `--megahit_fix_cpu_1`). See the [parameters page](https://nf-co.re/mag/dev/parameters).

The table is ordered first by 'memory', then 'Source' of defaults by specificity of definition, and then 'Module name' by alphabetical order.

| Module Name                          | CPU | Memory | Time | Source                |
| ------------------------------------ | --- | ------ | ---- | --------------------- |
| GTDBTK_CLASSIFYWF                    | 10  | 140.GB | 12.h | Named customisation   |
| CATPACK_BINS                         | 6   | 120.GB | 8.h  | Named customisation   |
| FLYE                                 | 12  | 72.GB  | 24.h | Named customisation   |
| METAMDBG_ASM                         | 12  | 72.GB  | 24.h | Named customisation   |
| FILTLONG                             | 8   | 64.GB  | 24.h | Named customisation   |
| CATPACK_CONTIGS                      | 6   | 60.GB  | 8.h  | Named customisation   |
| CHECKM_LINEAGEWF                     | 6   | 42.GB  | 8.h  | Named customisation   |
| MEGAHIT                              | 8   | 40.GB  | 16.h | Named customisation   |
| METASPADES                           | 10  | 40.GB  | 16.h | Named customisation   |
| METASPADESHYBRID                     | 10  | 40.GB  | 16.h | Named customisation   |
| PORECHOP_PORECHOP                    | 4   | 30.GB  | 4.h  | Named customisation   |
| BOWTIE2_HOST_REMOVAL_BUILD           | 10  | 20.GB  | 4.h  | Named customisation   |
| METABAT2_METABAT2                    | 8   | 20.GB  | 8.h  | Named customisation   |
| MAG_DEPTHS                           | 1   | 16.GB  | 4.h  | Named customisation   |
| BUSCO_BUSCO                          | 10  | 12.GB  | 8.h  | Named customisation   |
| BOWTIE2_HOST_REMOVAL_ALIGN           | 10  | 10.GB  | 6.h  | Named customisation   |
| NANOLYSE                             | 2   | 10.GB  | 3.h  | Named customisation   |
| BOWTIE2_ASSEMBLY_ALIGN               | 2   | 8.GB   | 8.h  | Named customisation   |
| BOWTIE2_PHIX_REMOVAL_ALIGN           | 4   | 8.GB   | 6.h  | Named customisation   |
| BOWTIE2_ASSEMBLY_BUILD               | 1   | 7.GB   | 4.h  | Named customisation   |
| CONCAT_BUSCO_TSV                     | 1   | 1.GB   | 4.h  | Named customisation   |
| CONCAT_CHECKM_TSV                    | 1   | 1.GB   | 4.h  | Named customisation   |
| CONCAT_CHECKM2_TSV                   | 1   | 1.GB   | 4.h  | Named customisation   |
| CONCAT_GUNC_CHECKM_TSV               | 1   | 1.GB   | 4.h  | Named customisation   |
| CONCAT_GUNC_TSV                      | 1   | 1.GB   | 4.h  | Named customisation   |
| COMEBIN_RUNCOMEBIN                   | 12  | 72.GB  | 16.h | Label: process_high   |
| CONCOCT_CONCOCT                      | 12  | 72.GB  | 16.h | Label: process_high   |
| GENOMAD_ENDTOEND                     | 12  | 72.GB  | 16.h | Label: process_high   |
| MINIMAP2_ALIGN                       | 12  | 72.GB  | 16.h | Label: process_high   |
| ADAPTERREMOVAL                       | 6   | 36.GB  | 8.h  | Label: process_medium |
| BBMAP_BBNORM                         | 6   | 36.GB  | 8.h  | Label: process_medium |
| BCFTOOLS_CONSENSUS                   | 6   | 36.GB  | 8.h  | Label: process_medium |
| BCFTOOLS_VIEW                        | 6   | 36.GB  | 8.h  | Label: process_medium |
| CATPACK_PREPARE                      | 6   | 36.GB  | 8.h  | Label: process_medium |
| CHECKM2_PREDICT                      | 6   | 36.GB  | 8.h  | Label: process_medium |
| CHOPPER                              | 6   | 36.GB  | 8.h  | Label: process_medium |
| DASTOOL_DASTOOL                      | 6   | 36.GB  | 8.h  | Label: process_medium |
| FASTP                                | 6   | 36.GB  | 8.h  | Label: process_medium |
| GUNC_RUN                             | 6   | 36.GB  | 8.h  | Label: process_medium |
| MAXBIN2                              | 6   | 36.GB  | 8.h  | Label: process_medium |
| METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS | 6   | 36.GB  | 8.h  | Label: process_medium |
| METABINNER_METABINNER                | 6   | 36.GB  | 8.h  | Label: process_medium |
| METAEUK_EASYPREDICT                  | 6   | 36.GB  | 8.h  | Label: process_medium |
| MMSEQS_DATABASES                     | 6   | 36.GB  | 8.h  | Label: process_medium |
| PORECHOP_ABI                         | 6   | 36.GB  | 8.h  | Label: process_medium |
| PYDAMAGE_ANALYZE                     | 6   | 36.GB  | 8.h  | Label: process_medium |
| SEMIBIN_SINGLEEASYBIN                | 6   | 36.GB  | 8.h  | Label: process_medium |
| TIARA_TIARA                          | 6   | 36.GB  | 8.h  | Label: process_medium |
| TRIMMOMATIC                          | 6   | 36.GB  | 8.h  | Label: process_medium |
| ADJUST_MAXBIN2_EXT                   | 2   | 12.GB  | 4.h  | Label: process_low    |
| BCFTOOLS_INDEX                       | 2   | 12.GB  | 4.h  | Label: process_low    |
| CHECKM_QA                            | 2   | 12.GB  | 4.h  | Label: process_low    |
| FASTQC                               | 2   | 12.GB  | 4.h  | Label: process_low    |
| FIND_CONCATENATE                     | 2   | 12.GB  | 4.h  | Label: process_low    |
| METABINNER_BINS                      | 2   | 12.GB  | 4.h  | Label: process_low    |
| METABINNER_KMER                      | 2   | 12.GB  | 4.h  | Label: process_low    |
| METABINNER_TOOSHORT                  | 2   | 12.GB  | 4.h  | Label: process_low    |
| MINIMAP2_INDEX                       | 2   | 12.GB  | 4.h  | Label: process_low    |
| NANOPLOT                             | 2   | 12.GB  | 4.h  | Label: process_low    |
| NANOQ                                | 2   | 12.GB  | 4.h  | Label: process_low    |
| PROKKA                               | 2   | 12.GB  | 4.h  | Label: process_low    |
| PYPOLCA_RUN                          | 2   | 12.GB  | 4.h  | Label: process_low    |
| RENAME_POSTDASTOOL                   | 2   | 12.GB  | 4.h  | Label: process_low    |
| RENAME_PREDASTOOL                    | 2   | 12.GB  | 4.h  | Label: process_low    |
| SAMTOOLS_INDEX                       | 2   | 12.GB  | 4.h  | Label: process_low    |
| SAMTOOLS_UNMAPPED                    | 2   | 12.GB  | 4.h  | Label: process_low    |
| SEQKIT_STATS                         | 2   | 12.GB  | 4.h  | Label: process_low    |
| SPLIT_FASTA                          | 2   | 12.GB  | 4.h  | Label: process_low    |
| SUMMARISE_PYDAMAGEBINS               | 2   | 12.GB  | 4.h  | Label: process_low    |
| ALE                                  | 1   | 6.GB   | 4.h  | Label: process_single |
| CAT_FASTQ                            | 1   | 6.GB   | 4.h  | Label: process_single |
| CATPACK_ADDNAMES                     | 1   | 6.GB   | 4.h  | Label: process_single |
| CATPACK_DOWNLOAD                     | 1   | 6.GB   | 4.h  | Label: process_single |
| CATPACK_SUMMARISE                    | 1   | 6.GB   | 4.h  | Label: process_single |
| CHECKM2_DATABASEDOWNLOAD             | 1   | 6.GB   | 4.h  | Label: process_single |
| CONCOCT_CONCOCTCOVERAGETABLE         | 1   | 6.GB   | 4.h  | Label: process_single |
| CONCOCT_CUTUPFASTA                   | 1   | 6.GB   | 4.h  | Label: process_single |
| CONCOCT_EXTRACTFASTABINS             | 1   | 6.GB   | 4.h  | Label: process_single |
| CONCOCT_MERGECUTUPCLUSTERING         | 1   | 6.GB   | 4.h  | Label: process_single |
| DASTOOL_FASTATOCONTIG2BIN            | 1   | 6.GB   | 4.h  | Label: process_single |
| FREEBAYES                            | 1   | 6.GB   | 4.h  | Label: process_single |
| GENOMAD_DOWNLOAD                     | 1   | 6.GB   | 4.h  | Label: process_single |
| GUNC_DOWNLOADDB                      | 1   | 6.GB   | 4.h  | Label: process_single |
| GUNC_MERGECHECKM                     | 1   | 6.GB   | 4.h  | Label: process_single |
| GUNZIP                               | 1   | 6.GB   | 4.h  | Label: process_single |
| MULTIQC                              | 1   | 6.GB   | 4.h  | Label: process_single |
| PRODIGAL                             | 1   | 6.GB   | 4.h  | Label: process_single |
| PYDAMAGE_FILTER                      | 1   | 6.GB   | 4.h  | Label: process_single |
| QSV_CAT                              | 1   | 6.GB   | 4.h  | Label: process_single |
| SAMTOOLS_FAIDX                       | 1   | 6.GB   | 4.h  | Label: process_single |
| SAMTOOLS_STATS                       | 1   | 6.GB   | 4.h  | Label: process_single |
| SEQTK_MERGEPE                        | 1   | 6.GB   | 4.h  | Label: process_single |
| TIARA_CLASSIFY                       | 1   | 6.GB   | 4.h  | Label: process_single |
| UNTAR                                | 1   | 6.GB   | 4.h  | Label: process_single |
| BIN_SUMMARY                          | 1   | 7.GB   | 4.h  | Default               |
| BOWTIE2_PHIX_REMOVAL_BUILD           | 1   | 7.GB   | 4.h  | Default               |
| CONVERT_DEPTHS                       | 1   | 7.GB   | 4.h  | Default               |
| GTDBTK_DB_PREPARATION                | 1   | 7.GB   | 4.h  | Default               |
| GTDBTK_SUMMARY                       | 1   | 7.GB   | 4.h  | Default               |
| MAG_DEPTHS_SUMMARY                   | 1   | 7.GB   | 4.h  | Default               |
| PREPARE_BIGMAG_SUMMARY               | 1   | 7.GB   | 4.h  | Default               |
| QUAST                                | 1   | 7.GB   | 4.h  | Default               |
| QUAST_BINS                           | 1   | 7.GB   | 4.h  | Default               |

_Table generated for nf-core/mag v5.5, using Claude Haiku 4.5, and corrected for accuracy by a human._

<!-- TODO mag: update when new modules added/defaults changed, if using an LLM, careful to tell it to follow module import aliases -->
