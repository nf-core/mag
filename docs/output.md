# nf-core/mag: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Quality control](#quality-control) of input reads - trimming and contaminant removal
* [Taxonomic classification of trimmed reads](#taxonomic-classification-of-trimmed-reads)
* [Assembly](#assembly) of trimmed reads
* [Binning](#binning) of assembled contigs
* [Taxonomic classification of binned genomes](#taxonomic-classification-of-binned-genomes)
* [Additional summary for binned genomes](#additional-summary-for-binned-genomes)
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Note that when specifying the parameter `--coassemble_group`, for the corresponding output filenames/directories of the assembly or downsteam processes the group ID, or more precisely the term `group-[group_id]`, will be used instead of the sample ID.

## Quality control

These steps trim away the adapter sequences present in input reads, trims away bad quality bases and sicard reads that are too short.
It also removes host contaminants and sequencing controls, such as PhiX or the Lambda phage.
FastQC is run for visualising the general quality metrics of the sequencing runs before and after trimming.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files:**

* `QC_shortreads/fastqc/`
  * `[sample]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files
  * `[sample].trimmed_[1/2]_fastqc.html`: FastQC report, containing quality metrics for trimmed and, if specified, filtered read files

### fastp

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the results folder and part of the MultiQC report.

**Output files:**

* `QC_shortreads/fastp/[sample]/`
  * `fastp.html`: Interactive report
  * `fastp.json`: Report in json format

### Remove PhiX sequences from short reads

The pipeline uses bowtie2 to map the reads against PhiX and removes mapped reads.

**Output files:**

* `QC_shortreads/remove_phix/`
  * `[sample].phix_removed.bowtie2.log`: Contains a brief log file indicating how many reads have been retained.

### Host read removal

The pipeline uses bowtie2 to map short reads against the host reference genome specified with `--host_genome` or `--host_fasta` and removes mapped reads. The information about discarded and retained reads is also included in the MultiQC report.

**Output files:**

* `QC_shortreads/remove_host/`
  * `[sample].host_removed.bowtie2.log`: Contains the bowtie2 log file indicating how many reads have been mapped as well as a file listing the read ids of discarded reads.

### Remove Phage Lambda sequences from long reads

The pipeline uses Nanolyse to map the reads against the Lambda phage and removes mapped reads.

**Output files:**

* `QC_longreads/NanoLyse/`
  * `[sample]_nanolyse.log`: Contains a brief log file indicating how many reads have been retained.

### Filtlong and porechop

The pipeline uses filtlong and porechop to perform quality control of the long reads that are eventually provided with the TSV input file.

No direct host read removal is performed for long reads.
However, since within this pipeline filtlong uses a read quality based on k-mer matches to the already filtered short reads, reads not overlapping those short reads might be discarded.
The lower the parameter `--longreads_length_weight`, the higher the impact of the read qualities for filtering.
For further documentation see the [filtlong online documentation](https://github.com/rrwick/Filtlong).

### Quality visualisation for long reads

NanoPlot is used to calculate various metrics and plots about the quality and length distribution of long reads. For more information about NanoPlot see the [online documentation](https://github.com/wdecoster/NanoPlot).

**Output files:**

* `QC_longreads/NanoPlot/[sample]/`
  * `raw_*.[png/html/txt]`: Plots and reports for raw data
  * `filtered_*.[png/html/txt]`: Plots and reports for filtered data

## Taxonomic classification of trimmed reads

### Kraken

Kraken2 classifies reads using a k-mer based approach as well as assigns taxonomy using a Lowest Common Ancestor (LCA) algorithm.

**Output files:**

* `Taxonomy/kraken2/[sample]/`
  * `kraken2.report`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
  * `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

### Centrifuge

Centrifuge is commonly used for the classification of DNA sequences from microbial samples. It uses an indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index.

More information on the [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) website

**Output files:**

* `Taxonomy/centrifuge/[sample]/`
  * `report.txt`: Tab-delimited result file. See the [centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output) for information about the fields
  * `kreport.txt`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
  * `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

## Assembly

Trimmed (short) reads are assembled with both megahit and SPAdes. Hybrid assembly is only supported by SPAdes.

### MEGAHIT

[MEGAHIT](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics short reads.

**Output files:**

* `Assembly/MEGAHIT/`
  * `[sample/group].contigs.fa.gz`: Compressed metagenome assembly in fasta format
  * `[sample/group].log`: Log file
  * `[sample/group]_QC/`: Directory containing QUAST files and Bowtie2 mapping logs
  * `[sample/group]_QC/MEGAHIT-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
  * `[sample/group]_QC/MEGAHIT-[sample/group]-[sampleToMap].bowtie2.log` and/or `[sample]_QC/MEGAHIT-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) was originally a single genome assembler that later added support for assembling metagenomes.

**Output files:**

* `Assembly/SPAdes/`
  * `[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
  * `[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
  * `[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
  * `[sample/group].log`: Log file
  * `[sample/group]_QC/`: Directory containing QUAST files and Bowtie2 mapping logs
  * `[sample/group]_QC/SPAdes-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
  * `[sample/group]_QC/SPAdes-[sample/group]-[sampleToMap].bowtie2.log` and/or `[sample]_QC/SPAdes-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

### SPAdesHybrid

SPAdesHybrid is a part of the [SPAdes](http://cab.spbu.ru/software/spades/) software and is used when the user provides both long and short reads.

**Output files:**

* `Assembly/SPAdesHybrid/`
  * `[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
  * `[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
  * `[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
  * `[sample/group].log`: Log file
  * `[sample/group]_QC/`: Directory containing QUAST files and Bowtie2 mapping logs
  * `[sample/group]_QC/SPAdesHybrid-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
  * `[sample/group]_QC/SPAdesHybrid-[sample/group]-[sampleToMap].bowtie2.log` and/or `[sample]_QC/SPAdesHybrid-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").

### Metagenome QC with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates metagenome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

**Output files:**

* `Assembly/[assembler]/[sample/group]_QC/`
  * `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
  * `quast.log`: QUAST log file
  * `predicted_genes/[assembler]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format

## Binning

### Contig sequencing depth

Sequencing depth per contig and sample is generated by `jgi_summarize_bam_contig_depths --outputDepth`. The values correspond to `(sum of exactely aligned bases) / ((contig length)-2*75)`. For example, for two reads aligned exactly with `10` and `9` bases on a 1000 bp long contig the depth is calculated by `(10+9)/(1000-2*75)` (1000bp length of contig minus 75bp from each end, which is excluded).

**Output files:**

* `GenomeBinning/`
  * `[assembler]-[sample/group]-depth.txt.gz`: Sequencing depth for each contig and sample or group, only for short reads.

### MetaBAT2

[MetaBAT2](https://bitbucket.org/berkeleylab/metabat) recovers genome bins (that is, contigs/scaffolds that all belongs to a same organism) from metagenome assemblies.

**Output files:**

* `GenomeBinning/MetaBAT2/`
  * `[assembler]-[sample/group].*.fa`: Genome bins retrieved from input assembly
  * `[assembler]-[sample/group].unbinned.*.fa`: Contigs that were not binned with other contigs but considered interesting. By default, these are at least 1 Mbp (`--min_length_unbinned_contigs`) in length and at most the 100 longest contigs (`--max_unbinned_contigs`) are reported

All the files and contigs in this folder will be assessed by QUAST and BUSCO.

**Output files:**

* `GenomeBinning/MetaBAT2/discarded/`
  * `*.lowDepth.fa.gz`: Low depth contigs that are filtered by MetaBat2
  * `*.tooShort.fa.gz`: Too short contigs that are filtered by MetaBat2
  * `*.unbinned.pooled.fa.gz`: Pooled unbinned contigs equal or above `--min_contig_size`, by default 1500 bp.
  * `*.unbinned.remaining.fa.gz`: Remaining unbinned contigs below `--min_contig_size`, by default 1500 bp, but not in any other file.

All the files in this folder contain small and/or unbinned contigs that are not further processed.

Files in these two folders contain all contigs of an assembly.

### Bin sequencing depth

For each genome bin the median sequencing depth is computed based on the corresponding contig depths given in `GenomeBinning/[assembler]-[sample/group]-depth.txt.gz`.

**Output files:**

* `GenomeBinning/`
  * `bin_depths_summary.tsv`: Summary of bin sequencing depths for all samples. Depths are available for samples mapped against the corresponding assembly, i.e. according to the mapping strategy specified with `--binning_map_mode`. Only for short reads.

### QC for metagenome assembled genomes with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates genome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

**Output files:**

* `GenomeBinning/QC/QUAST/[assembler]-[bin]/`
  * `report.*`: QUAST report in various formats, such as html, txt, tsv or tex
  * `quast.log`: QUAST log file
  * `predicted_genes/[assembler]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format
* `GenomeBinning/QC/`
  * `quast_summary.tsv`: QUAST output for all bins summarized

### QC for metagenome assembled genomes with BUSCO

[BUSCO](https://busco.ezlab.org/) is a tool used to assess the completeness of a genome assembly. It is run on all the genome bins and high quality contigs obtained by MetaBAT2. By default, BUSCO is run in automated lineage selection mode in which it first tries to select the domain and then a more specific lineage based on phylogenetic placement. If available, result files for both the selected domain lineage and the selected more specific lineage are placed in the output directory. If a lineage dataset is specified already with `--busco_reference`, only results for this specific lineage will be generated.

**Output files:**

* `GenomeBinning/QC/BUSCO/`
  * `[assembler]-[bin]_busco.log`: Log file containing the standard output of BUSCO.
  * `[assembler]-[bin]_busco.err`: File containing potential error messages returned from BUSCO.
  * `short_summary.domain.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results for the selected domain when run in automated lineage selection mode. Not available for bins for which a viral lineage was selected.
  * `short_summary.specific_lineage.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results in case a more specific lineage than the domain could be selected or for the lineage provided via `--busco_reference`.
  * `[assembler]-[bin]_buscos.[lineage].fna.gz`: Nucleotide sequence of all identified BUSCOs for used lineages (domain or specific).
  * `[assembler]-[bin]_buscos.[lineage].faa.gz`: Aminoacid sequence of all identified BUSCOs for used lineages (domain or specific).

If the parameter `--save_busco_reference` is set, additionally the used BUSCO lineage datasets are stored in the output directy.

**Output files:**

* `GenomeBinning/QC/BUSCO/`
  * `busco_downloads/`: All files and lineage datasets downloaded by BUSCO when run in automated lineage selection mode. (Can currently not be used to reproduce analysis, see the [nf-core/mag website documentation](https://nf-co.re/mag/usage#reproducibility) how to achieve reproducible BUSCO results).
  * `reference/*.tar.gz`: BUSCO reference lineage dataset that was provided via `--busco_reference`.

Besides the reference files or output files created by BUSCO, the following summary files will be generated:

**Output files:**

* `GenomeBinning/QC/`
  * `busco_summary.tsv`: A summary table of the BUSCO results, with % of marker genes found. If run in automated lineage selection mode, both the results for the selected domain and for the selected more specific lineage will be given, if available.

## Taxonomic classification of binned genomes

### CAT

[CAT](https://github.com/dutilh/CAT) is a toolkit for annotating contigs and bins from metagenome-assembled-genomes. The MAG pipeline uses CAT to assign taxonomy to genome bins based on the taxnomy of the contigs.

**Output files:**

* `Taxonomy/CAT/[assembler]/`
  * `[assembler]-[sample/group].ORF2LCA.names.txt.gz`: Tab-delimited files containing the lineage of each contig, with full lineage names
  * `[assembler]-[sample/group].bin2classification.names.txt.gz`: Taxonomy classification of the genome bins, with full lineage names
* `Taxonomy/CAT/[assembler]/raw/`
  * `[assembler]-[sample/group].concatenated.predicted_proteins.faa.gz`: Predicted protein sequences for each genome bin, in fasta format
  * `[assembler]-[sample/group].concatenated.predicted_proteins.gff.gz`: Predicted protein features for each genome bin, in gff format
  * `[assembler]-[sample/group].ORF2LCA.txt.gz`: Tab-delimited files containing the lineage of each contig
  * `[assembler]-[sample/group].bin2classification.txt.gz`: Taxonomy classification of the genome bins
  * `[assembler]-[sample/group].log`: Log files

If the parameters `--cat_db_generate` and `--save_cat_db` are set, additionally the generated CAT database is stored:

* `Taxonomy/CAT/CAT_prepare_*.tar.gz`: Generated and used CAT database.

### GTDB-Tk

[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) is a toolkit for assigning taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). nf-core/mag uses GTDB-Tk to classify binned genomes which satisfy certain quality criteria (i.e. completeness and contamination assessed with the BUSCO analysis).

**Output files:**

* `Taxonomy/GTDB-Tk/[assembler]-[sample/group]/`
  * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.summary.tsv`: Classifications for bacterial and archaeal genomes (see the [GTDB-Tk documentation for details](https://ecogenomics.github.io/GTDBTk/files/summary.tsv.html).
  * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.classify.tree.gz`: Reference tree in Newick format containing query genomes placed with pplacer.
  * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.markers_summary.tsv`: A summary of unique, duplicated, and missing markers within the 120 bacterial marker set, or the 122 archaeal marker set for each submitted genome.
  * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.msa.fasta.gz`: FASTA file containing MSA of submitted and reference genomes.
  * `gtdbtk.[assembler]-[sample/group].{bac120/ar122}.filtered.tsv`: A list of genomes with an insufficient number of amino acids in MSA.
  * `gtdbtk.[assembler]-[sample/group].*.log`: Log files.
  * `gtdbtk.[assembler]-[sample/group].failed_genomes.tsv`: A list of genomes for which the GTDB-Tk analysis failed, e.g. because Prodigal could not detect any genes.
* `Taxonomy/GTDB-Tk/gtdbtk_summary.tsv`: A summary table of the GTDB-Tk classification results for all bins, also containing bins which were discarded based on the BUSCO QC, which were filtered out by GTDB-Tk ((listed in `*.filtered.tsv`) or for which the analysis failed (listed in `*.failed_genomes.tsv`).

## Additional summary for binned genomes

**Output files:**

* `GenomeBinning/bin_summary.tsv`: Summary of bin sequencing depths together with BUSCO, QUAST and GTDB-Tk results, if at least one of the later was generated.

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
