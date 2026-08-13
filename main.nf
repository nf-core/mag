#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/mag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/mag
    Website: https://nf-co.re/mag
    Slack  : https://nfcore.slack.com/channels/mag
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAG                     } from './workflows/mag'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mag_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mag_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_mag_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_MAG {

    take:
        raw_short_reads  // channel: samplesheet read in from --input
        raw_long_reads
        input_assemblies

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MAG (
        raw_short_reads,  // channel: samplesheet read in from --input
        raw_long_reads,
        input_assemblies,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = MAG.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // iGenomes reference map, populated by conf/igenomes.config
    genomes: Map?

    // CSV samplesheet file containing information about the samples in the experiment.
    input: Path?

    // Specifies that the input is single-end reads.
    single_end: Boolean

    // Additional input CSV samplesheet containing information about pre-computed assemblies. When set, assembly is skipped and the supplied assemblies are used for downstream analysis.
    assembly_input: Path?

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String?

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // Do not load the iGenomes reference config.
    igenomes_ignore: Boolean

    // The base path to the igenomes reference files
    igenomes_base: String?

    // Git commit id for Institutional configs.
    custom_config_version: String?

    // Base directory for Institutional configs.
    custom_config_base: String?

    // Institutional config name.
    config_profile_name: String?

    // Institutional config description.
    config_profile_description: String?

    // Institutional config contact information.
    config_profile_contact: String?

    // Institutional config URL link.
    config_profile_url: String?

    // Display version and exit.
    version: Boolean

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String?

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String?

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String?

    // Display the help message. Accepts a parameter name to show help for that parameter only.
    help: String?

    // Display the full detailed help message.
    help_full: Boolean

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean

    // Fix number of CPUs for MEGAHIT to 1. Not increased with retries.
    megahit_fix_cpu_1: Boolean

    // Fix number of CPUs used by SPAdes. Not increased with retries.
    spades_fix_cpus: Integer = -1

    // Fix number of CPUs used by SPAdes hybrid. Not increased with retries.
    spadeshybrid_fix_cpus: Integer = -1

    // RNG seed for MetaBAT2.
    metabat_rng_seed: Integer = 1

    // Specify which adapter clipping tool to use.
    clip_tool: String = 'fastp'

    // Specify to save the resulting clipped FASTQ files to --outdir.
    save_clipped_reads: Boolean

    // The minimum length of reads must have to be retained for downstream analysis.
    reads_minlength: Integer = 15

    // Minimum phred quality value of a base to be qualified in fastp.
    fastp_qualified_quality: Integer = 15

    // The mean quality requirement used for per read sliding window cutting by fastp.
    fastp_cut_mean_quality: Integer = 15

    // Save reads that fail fastp filtering in a separate file. Not used downstream.
    fastp_save_trimmed_fail: Boolean

    // Turn on detecting and trimming of poly-G tails
    fastp_trim_polyg: Boolean

    // The minimum base quality for low-quality base trimming by AdapterRemoval.
    adapterremoval_minquality: Integer = 2

    // Turn on quality trimming by consecutive stretch of low quality bases, rather than by window.
    adapterremoval_trim_quality_stretch: Boolean

    // Forward read adapter to be trimmed by AdapterRemoval.
    adapterremoval_adapter1: String = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'

    // Reverse read adapter to be trimmed by AdapterRemoval for paired end data.
    adapterremoval_adapter2: String = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

    // Name of iGenomes reference for host contamination removal.
    host_genome: String?

    // Fasta reference file for host contamination removal.
    host_fasta: Path?

    // Bowtie2 index directory corresponding to `--host_fasta` reference file for host contamination removal.
    host_fasta_bowtie2index: Path?

    // Use the `--very-sensitive` instead of the`--sensitive`setting for Bowtie 2 to map reads against the host genome.
    host_removal_verysensitive: Boolean

    // Save the read IDs of removed host reads.
    host_removal_save_ids: Boolean

    // Specify to save input FASTQ files with host reads removed to --outdir.
    save_hostremoved_reads: Boolean

    // Keep reads similar to the Illumina internal standard PhiX genome.
    keep_phix: Boolean

    // Genome reference used to remove Illumina PhiX contaminant reads.
    phix_reference: Path?

    // Skip read preprocessing using fastp or adapterremoval.
    skip_clipping: Boolean

    // Skip all default QC steps for short reads (adapter trimming, phiX removal).
    skip_shortread_qc: Boolean

    // Skip running FastQC on short reads both before and after QC.
    skip_fastqc: Boolean

    // Specify to save input FASTQ files with phiX reads removed to --outdir.
    save_phixremoved_reads: Boolean

    // Run BBnorm to normalize sequence depth.
    bbnorm: Boolean

    // Set BBnorm target maximum depth to this number.
    bbnorm_target: Integer = 100

    // Set BBnorm minimum depth to this number.
    bbnorm_min: Integer = 5

    // Save normalized read files to output directory.
    save_bbnorm_reads: Boolean

    // Skip removing adapter sequences from long reads.
    skip_adapter_trimming: Boolean

    // Skip filtering long reads.
    skip_longread_filtering: Boolean

    // Skip all default QC steps for long reads (adapter trimming, filtering, removal of lambda sequences).
    skip_longread_qc: Boolean

    // Discard any read which is shorter than this value.
    longreads_min_length: Integer = 1000

    // Discard any read which has a mean quality score lower than this value.
    longreads_min_quality: Integer?

    // Keep this percent of bases. Only used by filtlong.
    longreads_keep_percent: Integer = 90

    // The higher the more important is read length when choosing the best reads. Only used by filtlong.
    longreads_length_weight: Integer = 10

    // Keep reads similar to the ONT internal standard Escherichia virus Lambda genome.
    keep_lambda: Boolean

    // Genome reference used to remove ONT Lambda contaminant reads.
    lambda_reference: Path?

    // Specify to save input FASTQ files with lamba reads removed  to --outdir.
    save_lambdaremoved_reads: Boolean

    // Specify to save the resulting clipped FASTQ files to --outdir.
    save_porechop_reads: Boolean

    // Specify to save the resulting length filtered long read FASTQ files to --outdir.
    save_filtered_longreads: Boolean

    // Specify which long read adapter trimming tool to use.
    longread_adaptertrimming_tool: String = 'porechop_abi'

    // Specify which long read filtering tool to use.
    longread_filtering_tool: String = 'chopper'

    // Filter long reads against short reads when using filtlong.
    filtlong_filtering_by_shortreads: Boolean

    // Database for taxonomic classification of metagenome assembled genomes. Can be either a zipped file or a directory containing the extracted output of such.
    cat_db: Path?

    // Generate CAT database.
    cat_db_generate: Boolean

    // Save the CAT database generated when specified by `--cat_db_generate`.
    save_cat_db: Boolean

    // Allow unofficial lineages in CAT classification.
    cat_allow_unofficial_lineages: Boolean

    // Classify unbinned contigs with CAT (contig mode).
    cat_classify_unbinned: Boolean

    // Specify to turn off CAT marking in output files most probable hits (when multiple) with an asterix.
    cat_no_suggestive_asterisks: Boolean

    // Skip the running of GTDB, as well as the automatic download of the database
    skip_gtdbtk: Boolean

    // Specify the location of a GTDBTK database. Can be either an uncompressed directory or a `.tar.gz` archive. If not specified will be downloaded for you when GTDBTK or binning QC is not skipped.
    gtdb_db: String = 'https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz'

    // Min. bin completeness (in %) required to apply GTDB-tk classification.
    gtdbtk_min_completeness: Float = 50.0

    // Max. bin contamination (in %) allowed to apply GTDB-tk classification.
    gtdbtk_max_contamination: Float = 10.0

    // Min. fraction of AA (in %) in the MSA for bins to be kept.
    gtdbtk_min_perc_aa: Float = 10.0

    // Min. alignment fraction to consider closest genome.
    gtdbtk_min_af: Float = 0.65

    // Number of CPUs used for the by GTDB-Tk run tool pplacer.
    gtdbtk_pplacer_cpus: Integer = 1

    // Speed up pplacer step of GTDB-Tk by loading to memory.
    gtdbtk_pplacer_useram: Boolean

    // Specify to have GTDBTk to use the full bacterial tree rather than the split tree (requires more memory!)
    gtdbtk_use_full_tree: Boolean

    // Specify to disable fast classification of genomes by ANI using skani in GTDB-Tk.
    gtdbtk_place_species: Boolean

    // [DEPRECATED] Use `--gtdbtk_place_species` instead. Specify to disable fast classification of genomes by ANI using skani in GTDB-Tk.
    gtdbtk_skip_aniscreen: Boolean

    // Run GTDB-Tk classification for all bins in a single job, rather than one job per sample/assembler/binner group.
    gtdbtk_single_job: Boolean

    // Co-assemble samples within one group, instead of assembling each sample separately.
    coassemble_group: Boolean

    // Additional custom options for SPAdes and SPAdesHybrid. Do not specify `--meta` as this will be added for you!
    spades_options: String?

    // Specify whether to use contigs or scaffolds assembled by SPAdes
    spades_downstreaminput: String = 'scaffolds'

    // Additional custom options for MEGAHIT.
    megahit_options: String?

    // Skip Illumina-only SPAdes assembly.
    skip_spades: Boolean

    // Skip SPAdes hybrid assembly.
    skip_spadeshybrid: Boolean

    // Skip MEGAHIT assembly.
    skip_megahit: Boolean

    // Skip ALE
    skip_ale: Boolean

    // Skip DeepMAsED assembly error detection (skips both features and predict)
    skip_deepmased: Boolean = true

    // Gzip DeepMAsED feature tables output. Useful for large assemblies to reduce disk usage.
    deepmased_features_gzip: Boolean

    // Random seed for numpy in DeepMAsED predict. Set for reproducible results.
    deepmased_predict_seed: Integer = 12

    // Enable ALE per-base output. This output can be very large (tens of GB).
    ale_per_base_output: Boolean

    // Skip metaQUAST.
    skip_quast: Boolean

    // Skip MetaDBG assembly.
    skip_metamdbg: Boolean

    // Skip Flye assembly.
    skip_flye: Boolean

    // Run PyPOLCA polishing on long-read assemblies before binning.
    run_pypolca: Boolean

    // Skip Prodigal gene prediction
    skip_prodigal: Boolean

    // Turn on Prokka complicance mode for truncating contig names for NCBI/ENA compatibility.
    prokka_with_compliance: Boolean

    // Specify sequencing centre name required for Prokka's compliance mode.
    prokka_compliance_centre: String?

    // Skip Prokka genome annotation.
    skip_prokka: Boolean

    // Specify to skip CDS/product searching in Prokka runs
    prokka_fast_mode: Boolean

    // [DEPRECATED] This parameter no longer has any effect and will be removed in a future release. MetaEuk only runs when `--metaeuk_db` or `--metaeuk_mmseqs_db` is supplied.
    skip_metaeuk: Boolean

    // A string containing the name of one of the databases listed in the [mmseqs2 documentation](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases). This database will be downloaded and formatted for eukaryotic genome annotation. Incompatible with --metaeuk_db.
    metaeuk_mmseqs_db: String?

    // Path to either a local fasta file of protein sequences, or to a directory containing an MMseqs2-formatted database, for annotation of eukaryotic genomes.
    metaeuk_db: Path?

    // Save the downloaded mmseqs2 database specified in `--metaeuk_mmseqs_db`.
    save_mmseqs_db: Boolean

    // Run virus identification.
    run_virus_identification: Boolean

    // Database for virus classification with geNomad
    genomad_db: Path?

    // Minimum geNomad score for a sequence to be considered viral
    genomad_min_score: Float = 0.7

    // Number of groups that geNomad's MMSeqs2 databse should be split into (reduced memory requirements)
    genomad_splits: Integer = 1

    // Defines mapping strategy to compute co-abundances for binning, i.e. which samples will be mapped against the assembly.
    binning_map_mode: String = 'group'

    // Skip metagenome binning entirely
    skip_binning: Boolean

    // Skip MetaBAT2 Binning
    skip_metabat2: Boolean

    // Skip MaxBin2 Binning
    skip_maxbin2: Boolean

    // Skip CONCOCT Binning
    skip_concoct: Boolean

    // Skip COMEBin Binning
    skip_comebin: Boolean

    // Skip MetaBinner Binning
    skip_metabinner: Boolean

    // Dataset scale for MetaBinner
    bin_metabinner_scale: String = 'large'

    // Skip SemiBin2 Binning
    skip_semibin: Boolean

    // RNG seed for SemiBin2.
    semibin_rng_seed: Integer = 1

    // Pre-trained model for SemiBin2 for single sample assemblies
    semibin_environment: String = 'global'

    // Minimum contig size to be considered for binning and for bin quality check.
    min_contig_size: Integer = 1500

    // Minimal length of contigs that are not part of any bin but treated as individual genome.
    min_length_unbinned_contigs: Integer = 1000000

    // Maximal number of contigs that are not part of any bin but treated as individual genome.
    max_unbinned_contigs: Integer = 100

    // Specify the shortest length a bin should be to retain for downstream processing (in base pairs)
    bin_min_size: Integer = 0

    // Specify the longest length a bin should be to retain for downstream processing (in base pairs). By default no limit.
    bin_max_size: Integer?

    // Limit the number of concurrent SEQKIT_STATS jobs used for bin size calculation.
    bin_seqkit_stats_max_forks: Integer?

    // Specify length of sub-contigs cut up prior CONCOCT binning
    bin_concoct_chunksize: Integer = 10000

    // Specify the overlap between each sub-contig prior CONCOCT binning
    bin_concoct_overlap: Integer = 0

    // Specify to not append the last contig less than sub-contig length to the last correct length contig
    bin_concoct_donotconcatlast: Boolean

    // Specify alternative Bowtie2 settings for aligning reads back against the assembly.
    bowtie2_mode: String?

    // Save the output of mapping raw reads back to assembled contigs
    save_assembly_mapped_reads: Boolean

    // Enable domain-level (prokaryote or eukaryote) classification of bins using Tiara. Processes which are domain-specific will then only receive bins matching the domain requirement.
    bin_domain_classification: Boolean

    // Specify which tool to use for domain classification of bins. Currently only 'tiara' is implemented.
    bin_domain_classification_tool: String = 'tiara'

    // Minimum contig length for Tiara to use for domain classification. For accurate classification, should be longer than 3000 bp.
    tiara_min_length: Integer = 3000

    // Exclude unbinned contigs in the post-binning steps (bin QC, taxonomic classification, and annotation steps).
    exclude_unbins_from_postbinning: Boolean

    // Specify a minimum percent identity filter for long reads mapping back to assembled contigs.
    longread_percentidentity: Float?

    // Specify a minimum percent identity filter for short reads mapping back to assembled contigs.
    shortread_percentidentity: Float?

    // Disable bin QC with BUSCO, CheckM or CheckM2.
    skip_binqc: Boolean

    // Enable running BUSCO during bin QC.
    run_busco: Boolean = true

    // Enable running CheckM during bin QC.
    run_checkm: Boolean

    // Enable running CheckM2 during bin QC.
    run_checkm2: Boolean

    // Download URL, local tar.gz archive, or local uncompressed directory for an *_odb10 or *_odb12 BUSCO lineage dataset.
    busco_db: Path?

    // Name of the BUSCO *_odb10 or *_odb12 lineage to check against. Additionally supports 'auto', 'auto_prok' and 'auto_euk' for automatic lineage selection mode.
    busco_db_lineage: String = 'auto'

    // Save the used BUSCO lineage datasets provided via `--busco_db`.
    save_busco_db: Boolean

    // Enable clean-up of temporary files created during BUSCO runs.
    busco_clean: Boolean

    // URL pointing to checkM database for auto download, if local path not supplied.
    checkm_download_url: String = 'https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz'

    // Path to local folder containing already downloaded and uncompressed CheckM database.
    checkm_db: Path?

    // Save the used CheckM reference files downloaded when not using --checkm_db parameter.
    save_checkm_data: Boolean

    // Path to local file of an already downloaded and uncompressed CheckM2 database file (.dmnd file).
    checkm2_db: Path?

    // CheckM2 database version number to download (Zenodo record ID, for reference check the canonical reference https://zenodo.org/records/5571251, and pick the Zenodo ID of the database version of your choice).
    checkm2_db_version: Integer = 14897628

    // Save the used CheckM2 reference files downloaded when not using --checkm2_db parameter.
    save_checkm2_data: Boolean

    // Turn on bin refinement using DAS Tool.
    refine_bins_dastool: Boolean

    // Specify single-copy gene score threshold for bin refinement.
    refine_bins_dastool_threshold: Float = 0.5

    // Specify to save contig to bin maps used for bin refinement
    refine_bins_dastool_savecontig2bin: Boolean

    // Specify which binning output is sent for downstream annotation, taxonomic classification, bin quality control etc.
    postbinning_input: String = 'raw_bins_only'

    // Turn on GUNC genome chimerism checks
    run_gunc: Boolean

    // Specify a path to a pre-downloaded GUNC dmnd database file
    gunc_db: Path?

    // Specify which database to auto-download if not supplying own
    gunc_database_type: String = 'progenomes_2.1'

    // Save the used GUNC reference files downloaded when not using --gunc_db parameter.
    gunc_save_db: Boolean

    // Make a BIgMAG input file including GUNC results.
    generate_bigmag_file: Boolean

    // Turn on/off the ancient DNA subworkflow
    ancient_dna: Boolean

    // PyDamage accuracy threshold
    pydamage_accuracy: Float = 0.5

    // deactivate damage correction of ancient contigs using variant and consensus calling
    skip_ancient_damagecorrection: Boolean

    // Ploidy for variant calling
    freebayes_ploidy: Integer = 1

    // minimum base quality required for variant calling
    freebayes_min_basequality: Integer = 20

    // minimum minor allele frequency for considering variants
    freebayes_minallelefreq: Float = 0.33

    // minimum genotype quality for considering a variant high quality
    bcftools_view_high_variant_quality: Integer = 30

    // minimum genotype quality for considering a variant medium quality
    bcftools_view_medium_variant_quality: Integer = 20

    // minimum number of bases supporting the alternative allele
    bcftools_view_minimal_allelesupport: Integer = 3
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MAG (
        PIPELINE_INITIALISATION.out.raw_short_reads,
        PIPELINE_INITIALISATION.out.raw_long_reads,
        PIPELINE_INITIALISATION.out.input_assemblies
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_MAG.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
