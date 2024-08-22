# nf-core/mag: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Quality control](#quality-control) of input reads - trimming and contaminant removal
- [Taxonomic classification of trimmed reads](#taxonomic-classification-of-trimmed-reads)
- [Digital sequencing normalisation](#digital-normalization-with-BBnorm)
- [Assembly](#assembly) of trimmed reads
- [Protein-coding gene prediction](#gene-prediction) of assemblies
- [Virus identification](#virus-identification-in-assemblies) of assemblies
- [Binning and binning refinement](#binning-and-binning-refinement) of assembled contigs
- [Taxonomic classification of binned genomes](#taxonomic-classification-of-binned-genomes)
- [Genome annotation of binned genomes](#genome-annotation-of-binned-genomes)
- [Additional summary for binned genomes](#additional-summary-for-binned-genomes)
- [Ancient DNA](#ancient-dna)
- [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Note that when specifying the parameter `--coassemble_group`, for the corresponding output filenames/directories of the assembly or downsteam processes the group ID, or more precisely the term `group-[group_id]`, will be used instead of the sample ID.

## Quality control

These steps trim away the adapter sequences present in input reads, trims away bad quality bases and sicard reads that are too short.
It also removes host contaminants and sequencing controls, such as PhiX or the Lambda phage.
FastQC is run for visualising the general quality metrics of the sequencing runs before and after trimming.

<!-- TODO: Add example MultiQC plots generated for this pipeline -->

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `QC_shortreads/fastqc/`
  - `[sample]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files
  - `[sample].trimmed_[1/2]_fastqc.html`: FastQC report, containing quality metrics for trimmed and, if specified, filtered read files

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### fastp

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads. Its output is in the results folder and part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

- `QC_shortreads/fastp/[sample]/`
  - `fastp.html`: Interactive report
  - `fastp.json`: Report in json format

</details>

### AdapterRemoval2

[AdapterRemoval](https://adapterremoval.readthedocs.io/en/stable/) searches for and removes remnant adapter sequences from High-Throughput Sequencing (HTS) data and (optionally) trims low quality bases from the 3' end of reads following adapter removal. It is popular in the field of palaeogenomics. The output logs are stored in the results folder, and as a part of the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

- `QC_shortreads/adapterremoval/[sample]/`
  - `[sample]_ar2.settings`: AdapterRemoval log file.

</details>

### Remove PhiX sequences from short reads

The pipeline uses bowtie2 to map the reads against PhiX and removes mapped reads.

<details markdown="1">
<summary>Output files</summary>

- `QC_shortreads/remove_phix/`
  - `[sample].phix_removed.bowtie2.log`: Contains a brief log file indicating how many reads have been retained.

</details>

### Host read removal

The pipeline uses bowtie2 to map short reads against the host reference genome specified with `--host_genome` or `--host_fasta` and removes mapped reads. The information about discarded and retained reads is also included in the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

- `QC_shortreads/remove_host/`
  - `[sample].host_removed.bowtie2.log`: Contains the bowtie2 log file indicating how many reads have been mapped.
  - `[sample].host_removed.mapped*.read_ids.txt`: Contains a file listing the read ids of discarded reads.

</details>

### Remove Phage Lambda sequences from long reads

The pipeline uses Nanolyse to map the reads against the Lambda phage and removes mapped reads.

<details markdown="1">
<summary>Output files</summary>

- `QC_longreads/NanoLyse/`
  - `[sample]_[run]_lambdafiltered.nanolyse.log`: Contains a brief log file indicating how many reads have been removed.

</details>

### Filtlong and porechop

The pipeline uses filtlong and porechop to perform quality control of the long reads that are eventually provided with the TSV input file.

No direct host read removal is performed for long reads.
However, since within this pipeline filtlong uses a read quality based on k-mer matches to the already filtered short reads, reads not overlapping those short reads might be discarded.
The lower the parameter `--longreads_length_weight`, the higher the impact of the read qualities for filtering.
For further documentation see the [filtlong online documentation](https://github.com/rrwick/Filtlong).

### Quality visualisation for long reads

NanoPlot is used to calculate various metrics and plots about the quality and length distribution of long reads. For more information about NanoPlot see the [online documentation](https://github.com/wdecoster/NanoPlot).

<details markdown="1">
<summary>Output files</summary>

- `QC_longreads/NanoPlot/[sample]/`
  - `raw_*.[png/html/txt]`: Plots and reports for raw data
  - `filtered_*.[png/html/txt]`: Plots and reports for filtered data

</details>

## Digital normalization with BBnorm

If the pipeline is called with the `--bbnorm` option, it will normalize sequencing depth of libraries prior assembly by removing reads to 1) reduce coverage of very abundant kmers and 2) delete very rare kmers (see `--bbnorm_target` and `--bbnorm_min` parameters).
When called in conjunction with `--coassemble_group`, BBnorm will operate on interleaved (merged) FastQ files, producing only a single output file.
If the `--save_bbnorm_reads` parameter is set, the resulting FastQ files are saved together with log output.

<details markdown="1">
<summary>Output files</summary>

- `bbmap/bbnorm/[sample]\*.fastq.gz`
- `bbmap/bbnorm/log/[sample].bbnorm.log`

</details>

## Taxonomic classification of trimmed reads

### Kraken

Kraken2 classifies reads using a k-mer based approach as well as assigns taxonomy using a Lowest Common Ancestor (LCA) algorithm.

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/kraken2/[sample]/`
  - `kraken2.report`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
  - `taxonomy.krona.html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

</details>

### Centrifuge

Centrifuge is commonly used for the classification of DNA sequences from microbial samples. It uses an indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index.

More information on the [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) website

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/centrifuge/[sample]/`
  - `[sample].kreport.txt`: Classification in the Kraken report format. See the [kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) for more details
  - `[sample].report.txt`: Tab-delimited result file. See the [centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output) for information about the fields
  - `[sample].results.txt`: Per read taxonomic classification information. See the [centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output) for more details
  - `[sample].html`: Interactive pie chart produced by [KronaTools](https://github.com/marbl/Krona/wiki)

</details>

## Assembly

Trimmed (short) reads are assembled with both megahit and SPAdes. Hybrid assembly is only supported by SPAdes.

### MEGAHIT

[MEGAHIT](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics short reads.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/MEGAHIT/`
  - `[sample/group].contigs.fa.gz`: Compressed metagenome assembly in fasta format
  - `[sample/group].log`: Log file
  - `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
    - `MEGAHIT-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
    - `MEGAHIT-[sample/group]-[sampleToMap].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").
    - `MEGAHIT-[sample].[bam/bai]`: Optionally saved BAM file of the Bowtie2 mapping of reads against the assembly.

</details>

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) was originally a single genome assembler that later added support for assembling metagenomes.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SPAdes/`
  - `[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
  - `[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
  - `[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
  - `[sample/group].log`: Log file
  - `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
    - `SPAdes-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
    - `SPAdes-[sample/group]-[sampleToMap].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").
    - `SPAdes-[sample].[bam/bai]`: Optionally saved BAM file of the Bowtie2 mapping of reads against the assembly.

</details>

### SPAdesHybrid

SPAdesHybrid is a part of the [SPAdes](http://cab.spbu.ru/software/spades/) software and is used when the user provides both long and short reads.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SPAdesHybrid/`
  - `[sample/group]_scaffolds.fasta.gz`: Compressed assembled scaffolds in fasta format
  - `[sample/group]_graph.gfa.gz`: Compressed assembly graph in gfa format
  - `[sample/group]_contigs.fasta.gz`: Compressed assembled contigs in fasta format
  - `[sample/group].log`: Log file
  - `QC/[sample/group]/`: Directory containing QUAST files and Bowtie2 mapping logs
    - `SPAdesHybrid-[sample].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the sample that the metagenome was assembled from, only present if `--coassemble_group` is not set.
    - `SPAdesHybrid-[sample/group]-[sampleToMap].bowtie2.log`: Bowtie2 log file indicating how many reads have been mapped from the respective sample ("sampleToMap").
    - `SPAdesHybrid-[sample].[bam/bai]`: Optionally saved BAM file of the Bowtie2 mapping of reads against the assembly.

</details>

### Metagenome QC with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates metagenome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/[assembler]/QC/[sample/group]/QUAST/`
  - `report.*`: QUAST report in various formats, such as html, pdf, tex, tsv, or txt
  - `transposed_report.*`: QUAST report that has been transposed into wide format (tex, tsv, or txt)
  - `quast.log`: QUAST log file
  - `metaquast.log`: MetaQUAST log file
  - `icarus.html`: Icarus main menu with links to interactive viewers
  - `icarus_viewers/contig_size_viewer.html`: Diagram of contigs that are ordered from longest to shortest
  - `basic_stats/cumulative_plot.pdf`: Shows the growth of contig lengths (contigs are ordered from largest to shortest)
  - `basic_stats/GC_content_plot.pdf`: Shows the distribution of GC content in the contigs
  - `basic_stats/[assembler]-[sample/group]_GC_content_plot.pdf`: Histogram of the GC percentage for the contigs
  - `basic_stats/Nx_plot.pdf`: Plot of Nx values as x varies from 0 to 100%.
  - `predicted_genes/[assembler]-[sample/group].rna.gff`: Contig positions for rRNA genes in gff version 3 format
  - `predicted_genes/barrnap.log`: Barrnap log file (ribosomal RNA predictor)

</details>

## Gene prediction

Protein-coding genes are predicted for each assembly.

<details markdown="1">
<summary>Output files</summary>

- `Annotation/Prodigal/`
  - `[assembler]-[sample/group].gff.gz`: Gene Coordinates in GFF format
  - `[assembler]-[sample/group].faa.gz`: The protein translation file consists of all the proteins from all the sequences in multiple FASTA format.
  - `[assembler]-[sample/group].fna.gz`: Nucleotide sequences of the predicted proteins using the DNA alphabet, not mRNA (so you will see 'T' in the output and not 'U').
  - `[assembler]-[sample/group]_all.txt.gz`: Information about start positions of genes.

</details>

## Virus identification in assemblies

### geNomad

[geNomad](https://github.com/apcamargo/genomad) identifies viruses and plasmids in sequencing data (isolates, metagenomes, and metatranscriptomes)

<details markdown="1">
<summary>Output files</summary>

- `VirusIdentification/geNomad/[assembler]-[sample/group]*/`
  - `[assembler]-[sample/group]*_annotate`
    - `[assembler]-[sample/group]*_taxonomy.tsv`: Taxonomic assignment data
  - `[assembler]-[sample/group]*_aggregated_classification`
    - `[assembler]-[sample/group]*_aggregated_classification.tsv`: Sequence classification in tabular format
  - `[assembler]-[sample/group]*_find_proviruses`
    - `[assembler]-[sample/group]*_provirus.tsv`: Characteristics of proviruses identified by geNomad
  - `[assembler]-[sample/group]*_summary`
    - `[assembler]-[sample/group]*_virus_summary.tsv`: Virus classification summary file in tabular format
    - `[assembler]-[sample/group]*_plasmid_summary.tsv`: Plasmid classification summary file in tabular format
    - `[assembler]-[sample/group]*_viruses_genes.tsv`: Virus gene annotation data in tabular format
    - `[assembler]-[sample/group]*_plasmids_genes.tsv`: Plasmid gene annotation data in tabular format
    - `[assembler]-[sample/group]*_viruses.fna`: Virus nucleotide sequences in FASTA format
    - `[assembler]-[sample/group]*_plasmids.fna`: Plasmid nucleotide sequences in FASTA format
    - `[assembler]-[sample/group]*_viruses_proteins.faa`: Virus protein sequences in FASTA format
    - `[assembler]-[sample/group]*_plasmids_proteins.faa`: Plasmid protein sequences in FASTA format
  - `[assembler]-[sample/group]*.log`: Plain text log file detailing the steps executed by geNomad (annotate, find-proviruses, marker-classification, nn-classification, aggregated-classification and summary)

</details>

## Binning and binning refinement

### Contig sequencing depth

Sequencing depth per contig and sample is generated by MetaBAT2's `jgi_summarize_bam_contig_depths --outputDepth`. The values correspond to `(sum of exactly aligned bases) / ((contig length)-2*75)`. For example, for two reads aligned exactly with `10` and `9` bases on a 1000 bp long contig the depth is calculated by `(10+9)/(1000-2*75)` (1000bp length of contig minus 75bp from each end, which is excluded).

These depth files are used for downstream binning steps.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/depths/contigs/`
  - `[assembler]-[sample/group]-depth.txt.gz`: Sequencing depth for each contig and sample or group, only for short reads.

</details>

### MetaBAT2

[MetaBAT2](https://bitbucket.org/berkeleylab/metabat) recovers genome bins (that is, contigs/scaffolds that all belongs to a same organism) from metagenome assemblies.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/MetaBAT2/`
  - `bins/[assembler]-[binner]-[sample/group].*.fa.gz`: Genome bins retrieved from input assembly
  - `unbinned/[assembler]-[binner]-[sample/group].unbinned.[1-9]*.fa.gz`: Contigs that were not binned with other contigs but considered interesting. By default, these are at least 1 Mbp (`--min_length_unbinned_contigs`) in length and at most the 100 longest contigs (`--max_unbinned_contigs`) are reported

</details>

All the files and contigs in these folders will be assessed by QUAST and BUSCO.

All other files that were discarded by the tool, or from the low-quality unbinned contigs, can be found here.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/MetaBAT2/discarded/`
  - `*.lowDepth.fa.gz`: Low depth contigs that are filtered by MetaBAT2
  - `*.tooShort.fa.gz`: Too short contigs that are filtered by MetaBAT2
- `GenomeBinning/MetaBAT2/unbinned/discarded/`
  - `*.unbinned.pooled.fa.gz`: Pooled unbinned contigs equal or above `--min_contig_size`, by default 1500 bp.
  - `*.unbinned.remaining.fa.gz`: Remaining unbinned contigs below `--min_contig_size`, by default 1500 bp, but not in any other file.

</details>

All the files in this folder contain small and/or unbinned contigs that are not further processed.

Files in these two folders contain all contigs of an assembly.

### MaxBin2

[MaxBin2](https://sourceforge.net/projects/maxbin2/) recovers genome bins (that is, contigs/scaffolds that all belongs to a same organism) from metagenome assemblies.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/MaxBin2/`
  - `bins/[assembler]-[binner]-[sample/group].*.fa.gz`: Genome bins retrieved from input assembly
  - `unbinned/[assembler]-[binner]-[sample/group].noclass.[1-9]*.fa.gz`: Contigs that were not binned with other contigs but considered interesting. By default, these are at least 1 Mbp (`--min_length_unbinned_contigs`) in length and at most the 100 longest contigs (`--max_unbinned_contigs`) are reported.

</details>

All the files and contigs in these folders will be assessed by QUAST and BUSCO.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/MaxBin2/discarded/`
  - `*.tooshort.gz`: Too short contigs that are filtered by MaxBin2
- `GenomeBinning/MaxBin2/unbinned/discarded/`
  - `*.noclass.pooled.fa.gz`: Pooled unbinned contigs equal or above `--min_contig_size`, by default 1500 bp.
  - `*.noclass.remaining.fa.gz`: Remaining unbinned contigs below `--min_contig_size`, by default 1500 bp, but not in any other file.

</details>

All the files in this folder contain small and/or unbinned contigs that are not further processed.

Files in these two folders contain all contigs of an assembly.

### CONCOCT

[CONCOCT](https://github.com/BinPro/CONCOCT) performs unsupervised binning of metagenomic contigs by using nucleotide composition, coverage data in multiple samples and linkage data from paired end reads.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/CONCOCT/`
  - `bins/[assembler]-[binner]-[sample/group].*.fa.gz`: Genome bins retrieved from input assembly
  - `stats/[assembler]-[binner]-[sample/group].csv`: Table indicating which contig goes with which cluster bin.
  - `stats/[assembler]-[binner]-[sample/group]*_gt1000.csv`: Various intermediate PCA statistics used for clustering.
  - `stats/[assembler]-[binner]-[sample/group]_*.tsv`: Coverage statistics of each sub-contig cut up by CONOCOCT prior in an intermediate step prior to binning. Likely not useful in most cases.
  - `stats/[assembler]-[binner]-[sample/group].log.txt`: CONCOCT execution log file.
  - `stats/[assembler]-[binner]-[sample/group]_*.args`: List of arguments used in CONCOCT execution.
  - </details>

All the files and contigs in these folders will be assessed by QUAST and BUSCO, if the parameter `--postbinning_input` is not set to `refined_bins_only`.

Note that CONCOCT does not output what it considers 'unbinned' contigs, therefore no 'discarded' contigs are produced here. You may still need to do your own manual curation of the resulting bins.

### DAS Tool

[DAS Tool](https://github.com/cmks/DAS_Tool) is an automated binning refinement method that integrates the results of a flexible number of binning algorithms to calculate an optimized, non-redundant set of bins from a single assembly. nf-core/mag uses this tool to attempt to further improve bins based on combining the MetaBAT2 and MaxBin2 binning output, assuming sufficient quality is met for those bins.

DAS Tool will remove contigs from bins that do not pass additional filtering criteria, and will discard redundant lower-quality output from binners that represent the same estimated 'organism', until the single highest quality bin is represented.

> ⚠️ If DAS Tool does not find any bins passing your selected threshold it will exit with an error. Such an error is 'ignored' by nf-core/mag, therefore you will not find files in the `GenomeBinning/DASTool/` results directory for that particular sample.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/DASTool/`
  - `[assembler]-[sample/group]_allBins.eval`: Tab-delimited description with quality and completeness metrics for the input bin sets. Quality and completeness are estimated by DAS TOOL using a scoring function based on the frequency of bacterial or archaeal reference single-copy genes (SCG). Please see note at the bottom of this section on file names.
  - `[assembler]-[sample/group]_DASTool_summary.tsv`: Tab-delimited description with quality and completeness metrics for the refined output bin sets.
  - `[assembler]-[sample/group]_DASTool_contig2bin.tsv`: File describing which contig is associated to which bin from the input binners.
  - `[assembler]-[sample/group]_DASTool.log`: Log file from the DAS Tool run describing the command executed and additional runtime information.
  - `[assembler]-[sample/group].seqlength`: Tab-delimited file describing the length of each contig.
  - `bins/[assembler]-[binner]Refined-[sample/group].*.fa`: Refined bins in fasta format.
  - `unbinned/[assembler]-DASToolUnbinned-[sample/group].*.fa`: Unbinned contigs from bin refinement in fasta format.

</details>

By default, only the raw bins (and unbinned contigs) from the actual binning methods, but not from the binning refinement with DAS Tool, will be used for downstream bin quality control, annotation and taxonomic classification. The parameter `--postbinning_input` can be used to change this behaviour.

⚠️ Due to ability to perform downstream QC of both raw and refined bins in parallel (via `--postbinning_input)`, bin names in DAS Tools's `*_allBins.eval` file will include `Refined`. However for this particular file, they _actually_ refer to the 'raw' input bins. The pipeline renames the input files prior to running DASTool to ensure they can be disambiguated from the original bin files in the downstream QC steps.

### Tiara

Tiara is a contig classifier that identifies the domain (prokarya, eukarya) of contigs within an assembly. This is used in this pipeline to rapidly and with few resources identify the most likely domain classification of each bin or unbin based on its contig identities.

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/Tiara/`
  - `[assembler]-[sample/group].tiara.txt` - Tiara output classifications (with probabilities) for all contigs within the specified sample/group assembly
  - `log_[assembler]-[sample/group].txt` - log file detailing the parameters used by the Tiara model for contig classification.
- `GenomeBinning/tiara_summary.tsv` - Summary of Tiara domain classification for all bins.

</details>

Typically, you would use `tiara_summary.tsv` as the primary file to see which bins or unbins have been classified to which domains at a glance, whereas the files in `Taxonomy/Tiara` provide classifications for each contig.

### Bin sequencing depth

For each bin or refined bin the median sequencing depth is computed based on the corresponding contig depths.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/depths/bins/`
  - `bin_depths_summary.tsv`: Summary of bin sequencing depths for all samples. Depths are available for samples mapped against the corresponding assembly, i.e. according to the mapping strategy specified with `--binning_map_mode`. Only for short reads.
  - `bin_refined_depths_summary.tsv`: Summary of sequencing depths for refined bins for all samples, if refinement was performed. Depths are available for samples mapped against the corresponding assembly, i.e. according to the mapping strategy specified with `--binning_map_mode`. Only for short reads.
  - `[assembler]-[binner]-[sample/group]-binDepths.heatmap.png`: Clustered heatmap showing bin abundances of the assembly across samples. Bin depths are transformed to centered log-ratios and bins as well as samples are clustered by Euclidean distance. Again, sample depths are available according to the mapping strategy specified with `--binning_map_mode`. If a sample produces only a single bin, a heatmap will not be provided.

</details>

### QC for metagenome assembled genomes with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates genome assemblies by computing various metrics. The QUAST output is in the bin directories shown below. This QUAST output is not shown in the MultiQC report.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/QUAST/[assembler]-[bin]/`
  - `report.*`: QUAST report in various formats, such as html, pdf, tex, tsv, or txt
  - `transposed_report.*`: QUAST report that has been transposed into wide format (tex, tsv, or txt)
  - `quast.log`: QUAST log file
  - `metaquast.log`: MetaQUAST log file
  - `icarus.html`: Icarus main menu with links to interactive viewers
  - `icarus_viewers/contig_size_viewer.html`: Diagram of contigs that are ordered from longest to shortest
  - `basic_stats/cumulative_plot.pdf`: Shows the growth of contig lengths (contigs are ordered from largest to shortest)
  - `basic_stats/GC_content_plot.pdf`: Shows the distribution of GC content in the contigs
  - `basic_stats/[assembler]-[bin]_GC_content_plot.pdf`: Histogram of the GC percentage for the contigs
  - `basic_stats/Nx_plot.pdf`: Plot of Nx values as x varies from 0 to 100%.
  - `predicted_genes/[assembler]-[bin].rna.gff`: Contig positions for rRNA genes in gff version 3 format
  - `predicted_genes/barrnap.log`: Barrnap log file (ribosomal RNA predictor)
- `GenomeBinning/QC/`
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group]-quast_summary.tsv`: QUAST output summarized per sample/condition.
  - `quast_summary.tsv`: QUAST output for all bins summarized

</details>

### QC for metagenome assembled genomes

#### BUSCO

[BUSCO](https://busco.ezlab.org/) is a tool used to assess the completeness of a genome assembly. It is run on all the genome bins and high quality contigs obtained by the applied binning and/or binning refinement methods (depending on the `--postbinning_input` parameter). By default, BUSCO is run in automated lineage selection mode in which it first tries to select the domain and then a more specific lineage based on phylogenetic placement. If available, result files for both the selected domain lineage and the selected more specific lineage are placed in the output directory. If a lineage dataset is specified already with `--busco_db`, only results for this specific lineage will be generated.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/BUSCO/`
  - `[assembler]-[bin]_busco.log`: Log file containing the standard output of BUSCO.
  - `[assembler]-[bin]_busco.err`: File containing potential error messages returned from BUSCO.
  - `short_summary.domain.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results for the selected domain when run in automated lineage selection mode. Not available for bins for which a viral lineage was selected.
  - `short_summary.specific_lineage.[lineage].[assembler]-[bin].txt`: BUSCO summary of the results in case a more specific lineage than the domain could be selected or for the lineage provided via `--busco_db`.
  - `[assembler]-[bin]_buscos.[lineage].fna.gz`: Nucleotide sequence of all identified BUSCOs for used lineages (domain or specific).
  - `[assembler]-[bin]_buscos.[lineage].faa.gz`: Aminoacid sequence of all identified BUSCOs for used lineages (domain or specific).
  - `[assembler]-[bin]_prodigal.gff`: Genes predicted with Prodigal.

</details>

If the parameter `--save_busco_db` is set, additionally the used BUSCO lineage datasets are stored in the output directory.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/BUSCO/`
  - `busco_downloads/`: All files and lineage datasets downloaded by BUSCO when run in automated lineage selection mode. (Can currently not be used to reproduce analysis, see the [nf-core/mag website documentation](https://nf-co.re/mag/usage#reproducibility) how to achieve reproducible BUSCO results).
  - `reference/*.tar.gz`: BUSCO reference lineage dataset that was provided via `--busco_db`.

</details>

Besides the reference files or output files created by BUSCO, the following summary files will be generated:

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/`
  - `busco_summary.tsv`: A summary table of the BUSCO results, with % of marker genes found. If run in automated lineage selection mode, both the results for the selected domain and for the selected more specific lineage will be given, if available.

</details>

#### CheckM

[CheckM](https://ecogenomics.github.io/CheckM/) CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes. It provides robust estimates of genome completeness and contamination by using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage

By default, nf-core/mag runs CheckM with the `check_lineage` workflow that places genome bins on a reference tree to define lineage-marker sets, to check for completeness and contamination based on lineage-specific marker genes. and then subsequently runs `qa` to generate the summary files.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/CheckM/`
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group]_qa.txt`: Detailed statistics about bins informing completeness and contamamination scores (output of `checkm qa`). This should normally be your main file to use to evaluate your results.
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group]_wf.tsv`: Overall summary file for completeness and contamination (output of `checkm lineage_wf`).
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group]/`: intermediate files for CheckM results, including CheckM generated annotations, log, lineage markers etc.
  - `checkm_summary.tsv`: A summary table of the CheckM results for all bins (output of `checkm qa`).

</details>

If the parameter `--save_checkm_reference` is set, additionally the used the CheckM reference datasets are stored in the output directory.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/CheckM/`
  - `checkm_downloads/`: All CheckM reference files downloaded from the CheckM FTP server, when not supplied by the user.
    - `checkm_data_2015_01_16/*`: a range of directories and files required for CheckM to run.

</details>

#### GUNC

[Genome UNClutterer (GUNC)](https://grp-bork.embl-community.io/gunc/index.html) is a tool for detection of chimerism and contamination in prokaryotic genomes resulting from mis-binning of genomic contigs from unrelated lineages. It does so by applying an entropy based score on taxonomic assignment and contig location of all genes in a genome. It is generally considered as a additional complement to CheckM results.

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/QC/gunc_summary.tsv`
- `GenomeBinning/QC/gunc_checkm_summary.tsv`
- `[gunc-database].dmnd`
- `GUNC/`
  - `raw/`
    - `[assembler]-[binner]-[domain]-[refinement]-[sample/group]/[fasta input file name]/GUNC_checkM.merged.tsv`: Per sample GUNC [output](https://grp-bork.embl-community.io/gunc/output.html) containing with taxonomic and completeness QC statistics.
  - `checkmmerged/`
    - `[assembler]-[binner]-[domain]-[refinement]-[sample/group]/[checkm input file name]/GUNC.progenomes_2.1.maxCSS_level.tsv`: Per sample GUNC output merged with output from [CheckM](#checkm)

</details>

GUNC will be run if specified with `--run_gunc` as a standalone, unless CheckM is also activated via `--qc_tool 'checkm'`, in which case GUNC output will be merged with the CheckM output using `gunc merge_checkm`.

If `--gunc_save_db` is specified, the output directory will also contain the requested database (progenomes, or GTDB) in DIAMOND format.

## Taxonomic classification of binned genomes

### CAT

[CAT](https://github.com/dutilh/CAT) is a toolkit for annotating contigs and bins from metagenome-assembled-genomes. The nf-core/mag pipeline uses CAT to assign taxonomy to genome bins based on the taxnomy of the contigs.

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/CAT/[assembler]/[binner]/`
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].ORF2LCA.names.txt.gz`: Tab-delimited files containing the lineage of each contig, with full lineage names
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].bin2classification.names.txt.gz`: Taxonomy classification of the genome bins, with full lineage names
- `Taxonomy/CAT/[assembler]/[binner]/raw/`
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].concatenated.predicted_proteins.faa.gz`: Predicted protein sequences for each genome bin, in fasta format
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].concatenated.predicted_proteins.gff.gz`: Predicted protein features for each genome bin, in gff format
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].ORF2LCA.txt.gz`: Tab-delimited files containing the lineage of each contig
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].bin2classification.txt.gz`: Taxonomy classification of the genome bins
  - `[assembler]-[binner]-[domain]-[refinement]-[sample/group].log`: Log files

</details>

If the parameters `--cat_db_generate` and `--save_cat_db` are set, additionally the generated CAT database is stored:

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/CAT/CAT_prepare_*.tar.gz`: Generated and used CAT database.

</details>

### GTDB-Tk

[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) is a toolkit for assigning taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). nf-core/mag uses GTDB-Tk to classify binned genomes which satisfy certain quality criteria (i.e. completeness and contamination assessed with the BUSCO analysis).

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/GTDB-Tk/[assembler]/[binner]/[sample/group]/`
  - `gtdbtk.[assembler]-[binner]-[sample/group].{bac120/ar122}.summary.tsv`: Classifications for bacterial and archaeal genomes (see the [GTDB-Tk documentation for details](https://ecogenomics.github.io/GTDBTk/files/summary.tsv.html)).
  - `gtdbtk.[assembler]-[binner]-[domain]-[refinement]-[sample/group].{bac120/ar122}.classify.tree.gz`: Reference tree in Newick format containing query genomes placed with pplacer.
  - `gtdbtk.[assembler]-[binner]-[domain]-[refinement]-[sample/group].{bac120/ar122}.markers_summary.tsv`: A summary of unique, duplicated, and missing markers within the 120 bacterial marker set, or the 122 archaeal marker set for each submitted genome.
  - `gtdbtk.[assembler]-[binner]-[domain]-[refinement]-[sample/group].{bac120/ar122}.msa.fasta.gz`: FASTA file containing MSA of submitted and reference genomes.
  - `gtdbtk.[assembler]-[binner]-[domain]-[refinement]-[sample/group].{bac120/ar122}.filtered.tsv`: A list of genomes with an insufficient number of amino acids in MSA.
  - `gtdbtk.[assembler]-[binner]-[domain]-[refinement]-[sample/group].*.log`: Log files.
  - `gtdbtk.[assembler]-[binner]-[domain]-[refinement]-[sample/group].failed_genomes.tsv`: A list of genomes for which the GTDB-Tk analysis failed, e.g. because Prodigal could not detect any genes.
- `Taxonomy/GTDB-Tk/gtdbtk_summary.tsv`: A summary table of the GTDB-Tk classification results for all bins, also containing bins which were discarded based on the BUSCO QC, which were filtered out by GTDB-Tk (listed in `*.filtered.tsv`) or for which the analysis failed (listed in `*.failed_genomes.tsv`).

</details>

## Genome annotation of binned genomes

### Prokka

Whole genome annotation is the process of identifying features of interest in a set of genomic DNA sequences, and labelling them with useful information. [Prokka](https://github.com/tseemann/prokka) is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files.

<details markdown="1">
<summary>Output files</summary>

- `Annotation/Prokka/[assembler]/[bin]/`
  - `[assembler]-[binner]-[bin].gff`: annotation in GFF3 format, containing both sequences and annotations
  - `[assembler]-[binner]-[bin].gbk`: annotation in GenBank format, containing both sequences and annotations
  - `[assembler]-[binner]-[bin].fna`: nucleotide FASTA file of the input contig sequences
  - `[assembler]-[binner]-[bin].faa`: protein FASTA file of the translated CDS sequences
  - `[assembler]-[binner]-[bin].ffn`: nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
  - `[assembler]-[binner]-[bin].sqn`: an ASN1 format "Sequin" file for submission to Genbank
  - `[assembler]-[binner]-[bin].fsa`: nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file
  - `[assembler]-[binner]-[bin].tbl`: feature Table file, used by "tbl2asn" to create the .sqn file
  - `[assembler]-[binner]-[bin].err`: unacceptable annotations - the NCBI discrepancy report.
  - `[assembler]-[binner]-[bin].log`: contains all the output that Prokka produced during its run
  - `[assembler]-[binner]-[bin].txt`: statistics relating to the annotated features found
  - `[assembler]-[binner]-[bin].tsv`: tab-separated file of all features (locus_tag, ftype, len_bp, gene, EC_number, COG, product)

</details>

### MetaEuk

In cases where eukaryotic genomes are recovered in binning, [MetaEuk](https://github.com/soedinglab/metaeuk) is also available to annotate eukaryotic genomes quickly with standards-compliant output files.

<details markdown="1">
<summary>Output files</summary>

- `Annotation/MetaEuk/[assembler]/[bin]`
  - `[assembler]-[binner]-[bin].fas`: fasta file of protein sequences identified by MetaEuk
  - `[assembler]-[binner]-[bin].codon.fas`: fasta file of nucleotide sequences corresponding to the protein sequences fasta
  - `[assembler]-[binner]-[bin].headersMap.tsv`: tab-separated table containing the information from each header in the fasta files
  - `[assembler]-[binner]-[bin].gff`: annotation in GFF3 format

</details>

## Additional summary for binned genomes

<details markdown="1">
<summary>Output files</summary>

- `GenomeBinning/bin_summary.tsv`: Summary of bin sequencing depths together with BUSCO, CheckM, QUAST and GTDB-Tk results, if at least one of the later was generated. This will also include refined bins if `--refine_bins_dastool` binning refinement is performed. Note that in contrast to the other tools, for CheckM the bin name given in the column "Bin Id" does not contain the ".fa" extension.

</details>

## Ancient DNA

Optional, only running when parameter `-profile ancient_dna` is specified.

### `PyDamage`

[Pydamage](https://github.com/maxibor/pydamage), is a tool to automate the process of ancient DNA damage identification and estimation from contigs. After modelling the ancient DNA damage using the C to T transitions, Pydamage uses a likelihood ratio test to discriminate between truly ancient, and modern contigs originating from sample contamination.

<details markdown="1">
<summary>Output files</summary>

- `Ancient_DNA/pydamage/analyze`
  - `[assembler]_[sample/group]/pydamage_results/pydamage_results.csv`: PyDamage raw result tabular file in `.csv` format. Format described here: [pydamage.readthedocs.io/en/0.62/output.html](https://pydamage.readthedocs.io/en/0.62/output.html)
- `Ancient_DNA/pydamage/filter`
  - `[assembler]_[sample/group]/pydamage_results/pydamage_results.csv`: PyDamage filtered result tabular file in `.csv` format. Format described here: [pydamage.readthedocs.io/en/0.62/output.html](https://pydamage.readthedocs.io/en/0.62/output.html)

</details>

### `variant_calling`

Because of aDNA damage, _de novo_ assemblers sometimes struggle to call a correct consensus on the contig sequence. To avoid this situation, the consensus is optionally re-called with a variant calling software using the reads aligned back to the contigs when `--run_ancient_damagecorrection` is supplied.

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/consensus`
  - `[assembler]_[sample/group].fa`: contigs sequence with re-called consensus from read-to-contig alignment
- `variant_calling/unfiltered`
  - `[assembler]_[sample/group].vcf.gz`: raw variant calls of the reads aligned back to the contigs.
- `variant_calling/filtered`
  - `[assembler]_[sample/group].filtered.vcf.gz`: quality filtered variant calls of the reads aligned back to the contigs.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

The general stats table at the top of the table will by default only display the most relevant pre- and post-processing statistics prior to assembly, i.e., FastQC, fastp/Adapter removal, and Bowtie2 PhiX and host removal mapping results.

Note that the FastQC raw and processed columns are right next to each other for improved visual comparability, however the processed columns represent the input reads _after_ fastp/Adapter Removal processing (the dedicated columns of which come directly after the two FastQC set of columns). Hover your cursor over each column name to see the which tool the column is derived from.

Summary tool-specific plots and tables of following tools are currently displayed (if activated):

- FastQC (pre- and post-trimming)
- fastp
- Adapter Removal
- bowtie2
- BUSCO
- QUAST
- Kraken2 / Centrifuge
- PROKKA

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
