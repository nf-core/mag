# Resource usage guidance

**nf-core/mag** is a bioinformatics best-practice analysis pipeline for the assembly, binning, and annotation of metagenomes.

Due to the nature of metagenomic samples, nf-core/mag can in some cases require large amounts of computational resources to execute.

In this page we provide information and guidance on how nf-core/mag works by default, and in some cases how to optimise computational resource usage.

## Grouping modules.

By default, nf-core/mag aims to make the runtime as efficiently as possible on HPC and similar infrastructure through generating more but shorter running jobs. 

For example, all pre-assembly steps typically run one job per input FASTQ file.
After assembly, each step runs on a per assembly or contig.
After binning, each step runs on a per-bin basis.

The only exception to this are report aggregation steps, which concatenate report tables from all samples together for a given step (e.g. CONCAT_QUAST_SUMMARY, or BIN_SUMMARY).
However these are computationally very light and should not impact a user.

## Default resource requests

The following table lists the default resource requests for every module in the pipeline.
Note that not all modules will run in every pipeline run, nor that the tool will necessarily actually use the requested amount.

_Table generated for nf-core/mag v5.5, using Claude Haiku 4.5, and verified for accuracy by a human_

| Module Name | CPU | Memory | Time | Source |
|-------------|-----|--------|------|--------|
| ADAPTERREMOVAL | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| ADJUST_MAXBIN2_EXT | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| ALE | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| BBMAP_BBNORM | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| BCFTOOLS_CONSENSUS | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| BCFTOOLS_INDEX | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| BCFTOOLS_VIEW | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| BIN_SUMMARY | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| BOWTIE2_ASSEMBLY_ALIGN | 2 * task.attempt | 8.GB * task.attempt | 8.h * task.attempt | Named process |
| BOWTIE2_ASSEMBLY_BUILD | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| BOWTIE2_REMOVAL_ALIGN | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| BOWTIE2_REMOVAL_BUILD | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| BUSCO_BUSCO | 10 * task.attempt | 12.GB * task.attempt | 8.h * task.attempt | Named process |
| CATPACK_ADDNAMES | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CATPACK_BINS | 6 * task.attempt | 120.GB * task.attempt | 8.h * task.attempt | Named process |
| CATPACK_CONTIGS | 6 * task.attempt | 60.GB * task.attempt | 8.h * task.attempt | Named process |
| CATPACK_DOWNLOAD | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CATPACK_PREPARE | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| CATPACK_SUMMARISE | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CAT_FASTQ | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CHECKM2_DATABASEDOWNLOAD | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CHECKM2_PREDICT | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| CHECKM_LINEAGEWF | 6 * task.attempt | 42.GB | 8.h * task.attempt | Named process |
| CHECKM_QA | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| CHOPPER | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| COMEBIN_RUNCOMEBIN | 12 * task.attempt | 72.GB * task.attempt | 16.h * task.attempt | Label: process_high |
| CONCOCT_CONCOCT | 12 * task.attempt | 72.GB * task.attempt | 16.h * task.attempt | Label: process_high |
| CONCOCT_CONCOCTCOVERAGETABLE | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CONCOCT_CUTUPFASTA | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CONCOCT_EXTRACTFASTABINS | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CONCOCT_MERGECUTUPCLUSTERING | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| CONVERT_DEPTHS | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| DASTOOL_DASTOOL | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| DASTOOL_FASTATOCONTIG2BIN | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| FASTP | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| FASTQC | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| FILTLONG | 8 * task.attempt | 64.GB * (2 ** (task.attempt - 1)) | 24.h * (2 ** (task.attempt - 1)) | Named process |
| FIND_CONCATENATE | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| FLYE | 12 * task.attempt | 72.GB * (2 ** (task.attempt - 1)) | 24.h * (2 ** (task.attempt - 1)) | Named process |
| FREEBAYES | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| GENOMAD_DOWNLOAD | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| GENOMAD_ENDTOEND | 12 * task.attempt | 72.GB * task.attempt | 16.h * task.attempt | Label: process_high |
| GTDBTK_CLASSIFYWF | 10 * task.attempt | 140.GB * task.attempt | 12.h * task.attempt | Named process |
| GTDBTK_DB_PREPARATION | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| GTDBTK_SUMMARY | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| GUNC_DOWNLOADDB | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| GUNC_MERGECHECKM | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| GUNC_RUN | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| GUNZIP | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| MAG_DEPTHS | 1 * task.attempt | 16.GB * task.attempt | 4.h * task.attempt | Named process |
| MAG_DEPTHS_SUMMARY | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| MAXBIN2 | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| MEGAHIT | params.megahit_fix_cpu_1 ? 1 : (8 * task.attempt) | 40.GB * task.attempt | 16.h * task.attempt | Named process |
| METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| METABAT2_METABAT2 | 8 * task.attempt | 20.GB * task.attempt | 8.h * task.attempt | Named process |
| METABINNER_BINS | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| METABINNER_KMER | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| METABINNER_METABINNER | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| METABINNER_TOOSHORT | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| METAEUK_EASYPREDICT | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| METAMDBG_ASM | 12 * task.attempt | 72.GB * (2 ** (task.attempt - 1)) | 24.h * (2 ** (task.attempt - 1)) | Named process |
| MINIMAP2_ALIGN | 12 * task.attempt | 72.GB * task.attempt | 16.h * task.attempt | Label: process_high |
| MINIMAP2_INDEX | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| MMSEQS_DATABASES | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| MULTIQC | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| NANOLYSE | 2 * task.attempt | 10.GB * task.attempt | 3.h * task.attempt | Named process |
| NANOPLOT | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| NANOQ | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| PORECHOP_ABI | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| PORECHOP_PORECHOP | 4 * task.attempt | 30.GB * task.attempt | 4.h * task.attempt | Named process |
| PREPARE_BIGMAG_SUMMARY | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| PRODIGAL | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| PROKKA | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| PYDAMAGE_ANALYZE | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| PYDAMAGE_FILTER | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| PYPOLCA_RUN | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| QSV_CAT | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| QUAST | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| QUAST_BINS | 1 * task.attempt | 7.GB * task.attempt | 4.h * task.attempt | Default |
| RENAME_POSTDASTOOL | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| RENAME_PREDASTOOL | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| SAMTOOLS_FAIDX | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| SAMTOOLS_INDEX | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| SAMTOOLS_STATS | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| SAMTOOLS_UNMAPPED | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| SEMIBIN_SINGLEEASYBIN | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| SEQKIT_STATS | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| SEQTK_MERGEPE | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| SPADES | 12 * task.attempt | 72.GB * task.attempt | 16.h * task.attempt | Label: process_high |
| SPLIT_FASTA | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| SUMMARISE_PYDAMAGEBINS | 2 * task.attempt | 12.GB * task.attempt | 4.h * task.attempt | Label: process_low |
| TIARA_CLASSIFY | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
| TIARA_TIARA | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| TRIMMOMATIC | 6 * task.attempt | 36.GB * task.attempt | 8.h * task.attempt | Label: process_medium |
| UNTAR | 1 | 6.GB * task.attempt | 4.h * task.attempt | Label: process_single |
