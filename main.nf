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
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // Input options
    input: Path?
    assembly_input: Path?
    single_end: Boolean

    // Short read preprocessing options
    skip_clipping: Boolean
    skip_shortread_qc: Boolean
    clip_tool: String = 'fastp'
    save_clipped_reads: Boolean
    reads_minlength: Integer = 15
    fastp_save_trimmed_fail: Boolean
    fastp_qualified_quality: Integer = 15
    fastp_cut_mean_quality: Integer = 15
    fastp_trim_polyg: Boolean
    adapterremoval_minquality: Integer = 2
    adapterremoval_trim_quality_stretch: Boolean
    adapterremoval_adapter1: String = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
    adapterremoval_adapter2: String = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    keep_phix: Boolean
    phix_reference: Path?
    save_phixremoved_reads: Boolean
    skip_fastqc: Boolean

    // Long read preprocessing options
    skip_longread_adapter_trimming: Boolean
    // [DEPRECATED] Use `--skip_longread_adapter_trimming` instead.
    skip_adapter_trimming: Boolean
    skip_longread_filtering: Boolean
    skip_longread_qc: Boolean
    longreads_min_length: Integer = 1000
    longreads_min_quality: Integer?
    longreads_keep_percent: Integer = 90
    longreads_length_weight: Integer = 10
    longread_adaptertrimming_tool: String = 'porechop_abi'
    longread_filtering_tool: String = 'chopper'
    filtlong_filtering_by_shortreads: Boolean
    save_porechop_reads: Boolean
    save_filtered_longreads: Boolean
    keep_lambda: Boolean
    lambda_reference: Path?
    save_lambdaremoved_reads: Boolean

    // Decontamination options
    host_fasta: Path?
    host_fasta_bowtie2index: Path?
    host_genome: String?
    host_removal_verysensitive: Boolean
    host_removal_save_ids: Boolean
    save_hostremoved_reads: Boolean

    // Other preprocessing options
    bbnorm: Boolean
    bbnorm_target: Integer = 100
    bbnorm_min: Integer = 5
    save_bbnorm_reads: Boolean

    // Assembly options
    skip_spades: Boolean
    spades_options: String?
    spades_downstreaminput: String = 'scaffolds'
    spades_fix_cpus: Integer = -1
    skip_spadeshybrid: Boolean
    spadeshybrid_fix_cpus: Integer = -1
    skip_megahit: Boolean
    megahit_options: String?
    megahit_fix_cpu_1: Boolean
    skip_flye: Boolean
    skip_metamdbg: Boolean
    coassemble_group: Boolean

    // Assembly QC and polishing options
    run_ale: Boolean
    skip_ale: Boolean
    run_deepmased: Boolean
    skip_quast: Boolean
    ale_per_base_output: Boolean
    deepmased_features_gzip: Boolean
    deepmased_predict_seed: Integer = 12
    run_pypolca: Boolean

    // Bin read mapping options
    binning_map_mode: String = 'group'
    bowtie2_mode: String?
    save_assembly_mapped_reads: Boolean
    shortread_percentidentity: Float?
    longread_percentidentity: Float?

    // Binning options
    skip_binning: Boolean
    min_contig_size: Integer = 1500
    max_contig_size: Integer?
    bin_min_size: Integer = 0
    bin_max_size: Integer?
    skip_comebin: Boolean
    skip_concoct: Boolean
    bin_concoct_chunksize: Integer = 10000
    bin_concoct_overlap: Integer = 0
    bin_concoct_donotconcatlast: Boolean
    skip_maxbin2: Boolean
    skip_metabat2: Boolean
    metabat_rng_seed: Integer = 1
    skip_metabinner: Boolean
    bin_metabinner_scale: String = 'large'
    skip_semibin: Boolean
    semibin_rng_seed: Integer = 1
    semibin_environment: String = 'global'
    min_length_unbinned_contigs: Integer = 1000000
    max_unbinned_contigs: Integer = 100
    bin_seqkit_stats_max_forks: Integer?
    exclude_unbins_from_postbinning: Boolean
    postbinning_input: String = 'raw_bins_only'

    // Bin refinement options
    refine_bins_dastool: Boolean
    refine_bins_dastool_threshold: Float = 0.5
    refine_bins_dastool_savecontig2bin: Boolean

    // Bin QC options
    skip_binqc: Boolean
    run_busco: Boolean = true
    busco_db: Path?
    busco_db_lineage: String = 'auto'
    save_busco_db: Boolean
    busco_clean: Boolean
    run_checkm: Boolean
    checkm_download_url: String = 'https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz'
    checkm_db: Path?
    save_checkm_data: Boolean
    run_checkm2: Boolean
    checkm2_db: Path?
    checkm2_db_version: Integer = 14897628
    save_checkm2_data: Boolean
    run_gunc: Boolean
    gunc_db: Path?
    gunc_database_type: String = 'progenomes_2.1'
    gunc_save_db: Boolean

    // Virus identification options
    run_virus_identification: Boolean
    genomad_db: Path?
    genomad_min_score: Float = 0.7
    genomad_splits: Integer = 1

    // Annotation options
    skip_prodigal: Boolean
    skip_prokka: Boolean
    prokka_fast_mode: Boolean
    prokka_with_compliance: Boolean
    prokka_compliance_centre: String?
    skip_bakta: Boolean
    bakta_db: Path?
    bakta_db_downloadtype: String = 'full'
    bakta_with_compliance: Boolean
    bakta_fast_mode: Boolean
    // [DEPRECATED] This parameter no longer has any effect and will be removed in a future release. MetaEuk only runs when `--metaeuk_db` or `--metaeuk_mmseqs_db` is supplied.
    skip_metaeuk: Boolean
    metaeuk_mmseqs_db: String?
    metaeuk_db: Path?
    save_mmseqs_db: Boolean

    // Taxonomy options
    bin_domain_classification: Boolean
    bin_domain_classification_tool: String = 'tiara'
    tiara_min_length: Integer = 3000
    cat_db: Path?
    cat_db_generate: Boolean
    save_cat_db: Boolean
    cat_allow_unofficial_lineages: Boolean
    cat_classify_unbinned: Boolean
    cat_no_suggestive_asterisks: Boolean
    skip_gtdbtk: Boolean
    gtdb_db: String = 'https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz'
    gtdbtk_min_completeness: Float = 50.0
    gtdbtk_max_contamination: Float = 10.0
    gtdbtk_min_perc_aa: Float = 10.0
    gtdbtk_min_af: Float = 0.65
    gtdbtk_pplacer_cpus: Integer = 1
    gtdbtk_pplacer_useram: Boolean
    gtdbtk_use_full_tree: Boolean
    gtdbtk_place_species: Boolean
    // [DEPRECATED] Use `--gtdbtk_place_species` instead. Specify to disable fast classification of genomes by ANI using skani in GTDB-Tk.
    gtdbtk_skip_aniscreen: Boolean
    gtdbtk_single_job: Boolean

    // Ancient DNA assembly validation options
    ancient_dna: Boolean
    pydamage_accuracy: Float = 0.5
    skip_ancient_damagecorrection: Boolean
    freebayes_ploidy: Integer = 1
    freebayes_min_basequality: Integer = 20
    freebayes_minallelefreq: Float = 0.33
    bcftools_view_high_variant_quality: Integer = 30
    bcftools_view_medium_variant_quality: Integer = 20
    bcftools_view_minimal_allelesupport: Integer = 3

    // MultiQC and BigMAG options
    multiqc_config: Path?
    multiqc_title: String?
    multiqc_logo: Path?
    max_multiqc_email_size: String = '25.MB'
    multiqc_methods_description: Path?
    generate_bigmag_file: Boolean
    dereplicate: Boolean
    dereplicate_annotate_all: Boolean
    dereplicate_min_completeness: Float = 50
    dereplicate_max_contamination: Float = 10
    dereplicate_ani: Float = 95

    // Reference genome options
    igenomes_ignore: Boolean
    igenomes_base: String?
    genomes: Map?

    // Boilerplate options
    outdir: String?
    publish_dir_mode: String?
    email: String?
    email_on_fail: String?
    plaintext_email: Boolean
    monochrome_logs: Boolean
    help: String?
    help_full: Boolean
    show_hidden: Boolean
    version: Boolean
    pipelines_testdata_base_path: String?
    trace_report_suffix: String?

    // Config options
    config_profile_name: String?
    config_profile_description: String?
    custom_config_version: String?
    custom_config_base: String?
    config_profile_contact: String?
    config_profile_url: String?

    // Schema validation default options
    validate_params: Boolean = true
}

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
    raw_short_reads // channel: samplesheet read in from --input
    raw_long_reads
    input_assemblies

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MAG(
        raw_short_reads,
        raw_long_reads,
        input_assemblies,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    multiqc_report              = MAG.out.multiqc_report // channel: /path/to/multiqc_report.html
    shortread_fastqc            = MAG.out.shortread_fastqc
    shortread_clipping          = MAG.out.shortread_clipping
    shortread_hostremoval       = MAG.out.shortread_hostremoval
    shortread_phixremoval       = MAG.out.shortread_phixremoval
    shortread_bbnorm            = MAG.out.shortread_bbnorm
    longread_adaptertrimming    = MAG.out.longread_adaptertrimming
    longread_filtering          = MAG.out.longread_filtering
    longread_nanolyse           = MAG.out.longread_nanolyse
    longread_nanoplot           = MAG.out.longread_nanoplot
    longread_hostremoval_index  = MAG.out.longread_hostremoval_index
    longread_hostremoval_bam    = MAG.out.longread_hostremoval_bam
    longread_hostremoval_fastq  = MAG.out.longread_hostremoval_fastq
    longread_hostremoval_stats  = MAG.out.longread_hostremoval_stats
    assembly                    = MAG.out.assembly
    assembly_pypolca            = MAG.out.assembly_pypolca
    assembly_quast              = MAG.out.assembly_quast
    assembly_mapping            = MAG.out.assembly_mapping
    assembly_mapping_full       = MAG.out.assembly_mapping_full
    assembly_ale                = MAG.out.assembly_ale
    assembly_deepmased_features = MAG.out.assembly_deepmased_features
    assembly_deepmased          = MAG.out.assembly_deepmased
    virus_identification        = MAG.out.virus_identification
    adna_freebayes              = MAG.out.adna_freebayes
    adna_filtered               = MAG.out.adna_filtered
    adna_consensus              = MAG.out.adna_consensus
    adna_pydamage_analyze       = MAG.out.adna_pydamage_analyze
    adna_pydamage_filter        = MAG.out.adna_pydamage_filter
    binning_depths_contigs      = MAG.out.binning_depths_contigs
    binning_metabat2            = MAG.out.binning_metabat2
    binning_maxbin2             = MAG.out.binning_maxbin2
    binning_concoct_stats       = MAG.out.binning_concoct_stats
    binning_concoct_bins        = MAG.out.binning_concoct_bins
    binning_comebin             = MAG.out.binning_comebin
    binning_metabinner          = MAG.out.binning_metabinner
    binning_semibin             = MAG.out.binning_semibin
    binning_unbinned            = MAG.out.binning_unbinned
    binning_dastool_contig2bin  = MAG.out.binning_dastool_contig2bin
    binning_dastool             = MAG.out.binning_dastool
    binning_dastool_bins        = MAG.out.binning_dastool_bins
    binning_depths_bins         = MAG.out.binning_depths_bins
    binning_dereplication       = MAG.out.binning_dereplication
    binning_summary             = MAG.out.binning_summary
    binning_bigmag              = MAG.out.binning_bigmag
    binning_contig2bin          = MAG.out.binning_contig2bin
    binqc_busco                 = MAG.out.binqc_busco
    binqc_busco_db              = MAG.out.binqc_busco_db
    binqc_checkm_db             = MAG.out.binqc_checkm_db
    binqc_checkm                = MAG.out.binqc_checkm
    binqc_checkm2_db            = MAG.out.binqc_checkm2_db
    binqc_checkm2               = MAG.out.binqc_checkm2
    binqc_gunc_db               = MAG.out.binqc_gunc_db
    binqc_gunc_raw              = MAG.out.binqc_gunc_raw
    binqc_gunc_checkm           = MAG.out.binqc_gunc_checkm
    binqc_summaries             = MAG.out.binqc_summaries
    binqc_pydamage_bins         = MAG.out.binqc_pydamage_bins
    binqc_tiara                 = MAG.out.binqc_tiara
    taxonomy_tiara              = MAG.out.taxonomy_tiara
    taxonomy_cat_db             = MAG.out.taxonomy_cat_db
    taxonomy_cat_bins           = MAG.out.taxonomy_cat_bins
    taxonomy_cat_summary        = MAG.out.taxonomy_cat_summary
    taxonomy_cat_unbins         = MAG.out.taxonomy_cat_unbins
    taxonomy_gtdbtk             = MAG.out.taxonomy_gtdbtk
    taxonomy_gtdbtk_summary     = MAG.out.taxonomy_gtdbtk_summary
    annotation_prodigal         = MAG.out.annotation_prodigal
    annotation_prokka           = MAG.out.annotation_prokka
    annotation_bakta            = MAG.out.annotation_bakta
    annotation_mmseqs_db        = MAG.out.annotation_mmseqs_db
    annotation_metaeuk          = MAG.out.annotation_metaeuk
    multiqc                     = MAG.out.multiqc
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MAG(
        PIPELINE_INITIALISATION.out.raw_short_reads,
        PIPELINE_INITIALISATION.out.raw_long_reads,
        PIPELINE_INITIALISATION.out.input_assemblies,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_MAG.out.multiqc_report,
    )

    publish:
    shortread_fastqc            = NFCORE_MAG.out.shortread_fastqc
    shortread_clipping          = NFCORE_MAG.out.shortread_clipping
    shortread_hostremoval       = NFCORE_MAG.out.shortread_hostremoval
    shortread_phixremoval       = NFCORE_MAG.out.shortread_phixremoval
    shortread_bbnorm            = NFCORE_MAG.out.shortread_bbnorm
    longread_adaptertrimming    = NFCORE_MAG.out.longread_adaptertrimming
    longread_filtering          = NFCORE_MAG.out.longread_filtering
    longread_nanolyse           = NFCORE_MAG.out.longread_nanolyse
    longread_nanoplot           = NFCORE_MAG.out.longread_nanoplot
    longread_hostremoval_index  = NFCORE_MAG.out.longread_hostremoval_index
    longread_hostremoval_bam    = NFCORE_MAG.out.longread_hostremoval_bam
    longread_hostremoval_fastq  = NFCORE_MAG.out.longread_hostremoval_fastq
    longread_hostremoval_stats  = NFCORE_MAG.out.longread_hostremoval_stats
    assembly                    = NFCORE_MAG.out.assembly
    assembly_pypolca            = NFCORE_MAG.out.assembly_pypolca
    assembly_quast              = NFCORE_MAG.out.assembly_quast
    assembly_mapping            = NFCORE_MAG.out.assembly_mapping
    assembly_mapping_full       = NFCORE_MAG.out.assembly_mapping_full
    assembly_ale                = NFCORE_MAG.out.assembly_ale
    assembly_deepmased_features = NFCORE_MAG.out.assembly_deepmased_features
    assembly_deepmased          = NFCORE_MAG.out.assembly_deepmased
    virus_identification        = NFCORE_MAG.out.virus_identification
    adna_freebayes              = NFCORE_MAG.out.adna_freebayes
    adna_filtered               = NFCORE_MAG.out.adna_filtered
    adna_consensus              = NFCORE_MAG.out.adna_consensus
    adna_pydamage_analyze       = NFCORE_MAG.out.adna_pydamage_analyze
    adna_pydamage_filter        = NFCORE_MAG.out.adna_pydamage_filter
    binning_depths_contigs      = NFCORE_MAG.out.binning_depths_contigs
    binning_metabat2            = NFCORE_MAG.out.binning_metabat2
    binning_maxbin2             = NFCORE_MAG.out.binning_maxbin2
    binning_concoct_stats       = NFCORE_MAG.out.binning_concoct_stats
    binning_concoct_bins        = NFCORE_MAG.out.binning_concoct_bins
    binning_comebin             = NFCORE_MAG.out.binning_comebin
    binning_metabinner          = NFCORE_MAG.out.binning_metabinner
    binning_semibin             = NFCORE_MAG.out.binning_semibin
    binning_unbinned            = NFCORE_MAG.out.binning_unbinned
    binning_dastool_contig2bin  = NFCORE_MAG.out.binning_dastool_contig2bin
    binning_dastool             = NFCORE_MAG.out.binning_dastool
    binning_dastool_bins        = NFCORE_MAG.out.binning_dastool_bins
    binning_depths_bins         = NFCORE_MAG.out.binning_depths_bins
    binning_dereplication       = NFCORE_MAG.out.binning_dereplication
    binning_summary             = NFCORE_MAG.out.binning_summary
    binning_bigmag              = NFCORE_MAG.out.binning_bigmag
    binning_contig2bin          = NFCORE_MAG.out.binning_contig2bin
    binqc_busco                 = NFCORE_MAG.out.binqc_busco
    binqc_busco_db              = NFCORE_MAG.out.binqc_busco_db
    binqc_checkm_db             = NFCORE_MAG.out.binqc_checkm_db
    binqc_checkm                = NFCORE_MAG.out.binqc_checkm
    binqc_checkm2_db            = NFCORE_MAG.out.binqc_checkm2_db
    binqc_checkm2               = NFCORE_MAG.out.binqc_checkm2
    binqc_gunc_db               = NFCORE_MAG.out.binqc_gunc_db
    binqc_gunc_raw              = NFCORE_MAG.out.binqc_gunc_raw
    binqc_gunc_checkm           = NFCORE_MAG.out.binqc_gunc_checkm
    binqc_summaries             = NFCORE_MAG.out.binqc_summaries
    binqc_pydamage_bins         = NFCORE_MAG.out.binqc_pydamage_bins
    binqc_tiara                 = NFCORE_MAG.out.binqc_tiara
    taxonomy_tiara              = NFCORE_MAG.out.taxonomy_tiara
    taxonomy_cat_db             = NFCORE_MAG.out.taxonomy_cat_db
    taxonomy_cat_bins           = NFCORE_MAG.out.taxonomy_cat_bins
    taxonomy_cat_summary        = NFCORE_MAG.out.taxonomy_cat_summary
    taxonomy_cat_unbins         = NFCORE_MAG.out.taxonomy_cat_unbins
    taxonomy_gtdbtk             = NFCORE_MAG.out.taxonomy_gtdbtk
    taxonomy_gtdbtk_summary     = NFCORE_MAG.out.taxonomy_gtdbtk_summary
    annotation_prodigal         = NFCORE_MAG.out.annotation_prodigal
    annotation_prokka           = NFCORE_MAG.out.annotation_prokka
    annotation_bakta            = NFCORE_MAG.out.annotation_bakta
    annotation_mmseqs_db        = NFCORE_MAG.out.annotation_mmseqs_db
    annotation_metaeuk          = NFCORE_MAG.out.annotation_metaeuk
    multiqc                     = NFCORE_MAG.out.multiqc
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    shortread_fastqc {
        path 'QC_shortreads/fastqc'
    }
    shortread_clipping {
        path { meta, _files -> "QC_shortreads/${params.clip_tool}/${meta.id}" }
    }
    shortread_hostremoval {
        path 'QC_shortreads/remove_host'
    }
    shortread_phixremoval {
        path 'QC_shortreads/remove_phix'
    }
    shortread_bbnorm {
        path { _meta, fastq, log ->
            fastq >> 'bbmap/bbnorm/'
            log >> 'bbmap/bbnorm/logs/'
        }
    }
    longread_adaptertrimming {
        path 'QC_longreads/porechop'
    }
    longread_filtering {
        path "QC_longreads/${[filtlong: 'Filtlong', nanoq: 'Nanoq', chopper: 'Chopper'][params.longread_filtering_tool]}"
    }
    longread_nanolyse {
        path 'QC_longreads/NanoLyse'
    }
    longread_nanoplot {
        path { meta, _files -> "QC_longreads/NanoPlot/${meta.id}" }
    }
    longread_hostremoval_index {
        path 'QC_longreads/minimap2/index'
    }
    longread_hostremoval_bam {
        path 'QC_longreads/minimap2/align'
    }
    longread_hostremoval_fastq {
        path 'QC_longreads/samtools/fastq'
    }
    longread_hostremoval_stats {
        path 'QC_longreads/samtools/stats'
    }
    assembly {
        path { meta, _files -> "Assembly/${meta.assembler}" }
    }
    assembly_pypolca {
        path { meta, polished, vcf, report ->
            polished >> "Assembly/${meta.assembler}-pypolca/${meta.assembler}-${meta.id}/${polished.name.replace('_corrected.fasta', '_pypolca.fasta')}"
            vcf >> "Assembly/${meta.assembler}-pypolca/${meta.assembler}-${meta.id}/"
            report >> "Assembly/${meta.assembler}-pypolca/${meta.assembler}-${meta.id}/"
        }
    }
    assembly_quast {
        path { meta, _files -> "Assembly/${meta.assembler}/QC/${meta.id}" }
    }
    assembly_mapping {
        path { meta, _files -> "Assembly/${meta.assembler}/QC/${meta.id}" }
    }
    assembly_mapping_full {
        path { meta, _files -> "Assembly/${meta.assembler}/QC/${meta.id}/full_assembly" }
    }
    assembly_ale {
        path { meta, _files -> "Assembly/${meta.assembler}/QC/${meta.id}/ALE" }
    }
    assembly_deepmased_features {
        path { meta, _files -> "Assembly/${meta.assembler}/QC/${meta.id}/DeepMAsED/features" }
    }
    assembly_deepmased {
        path { meta, _files -> "Assembly/${meta.assembler}/QC/${meta.id}/DeepMAsED" }
    }
    virus_identification {
        path { meta, _files -> "VirusIdentification/geNomad/${meta.id}" }
    }
    adna_freebayes {
        path 'Ancient_DNA/variant_calling/freebayes'
    }
    adna_filtered {
        path 'Ancient_DNA/variant_calling/filtered'
    }
    adna_consensus {
        path 'Ancient_DNA/variant_calling/consensus'
    }
    adna_pydamage_analyze {
        path { meta, _files -> "Ancient_DNA/pydamage/analyze/${meta.assembler}-${meta.id}" }
    }
    adna_pydamage_filter {
        path { meta, _files -> "Ancient_DNA/pydamage/filter/${meta.assembler}-${meta.id}" }
    }
    binning_depths_contigs {
        path 'GenomeBinning/depths/contigs'
    }
    binning_metabat2 {
        path { _meta, bins, tooshort, lowdepth ->
            bins >> 'GenomeBinning/MetaBAT2/bins/'
            tooshort >> 'GenomeBinning/MetaBAT2/discarded/'
            lowdepth >> 'GenomeBinning/MetaBAT2/discarded/'
        }
    }
    binning_maxbin2 {
        path { _meta, bins, tooshort, summary, abundance ->
            bins >> 'GenomeBinning/MaxBin2/bins/'
            tooshort >> 'GenomeBinning/MaxBin2/discarded/'
            summary >> 'GenomeBinning/MaxBin2/'
            abundance >> 'GenomeBinning/MaxBin2/'
        }
    }
    binning_concoct_stats {
        path 'GenomeBinning/CONCOCT/stats'
    }
    binning_concoct_bins {
        path { _meta, bin ->
            bin >> "GenomeBinning/CONCOCT/bins/${bin.name}"
        }
    }
    binning_comebin {
        path { _meta, bins, tsv, log, embeddings, covembeddings ->
            bins >> 'GenomeBinning/COMEBin/bins/'
            tsv >> 'GenomeBinning/COMEBin/stats/'
            log >> 'GenomeBinning/COMEBin/stats/'
            embeddings >> 'GenomeBinning/COMEBin/stats/'
            covembeddings >> 'GenomeBinning/COMEBin/stats/'
        }
    }
    binning_metabinner {
        path { _meta, membership, log, bins, unbinned ->
            membership >> 'GenomeBinning/MetaBinner/stats/'
            log >> 'GenomeBinning/MetaBinner/stats/'
            bins >> 'GenomeBinning/MetaBinner/'
            unbinned >> 'GenomeBinning/MetaBinner/unbinned/'
        }
    }
    binning_semibin {
        path { _meta, bins, log ->
            bins >> 'GenomeBinning/SemiBin2/'
            log >> 'GenomeBinning/SemiBin2/log/'
        }
    }
    binning_unbinned {
        path { meta, unbinned, pooled ->
            unbinned >> "GenomeBinning/${meta.binner}/unbinned/"
            pooled >> "GenomeBinning/${meta.binner}/unbinned/discarded/"
        }
    }
    binning_dastool_contig2bin {
        path 'GenomeBinning/DASTool/fastatocontig2bin'
    }
    binning_dastool {
        path 'GenomeBinning/DASTool'
    }
    binning_dastool_bins {
        path { _meta, bins, unbins ->
            bins >> 'GenomeBinning/DASTool/bins/'
            unbins >> 'GenomeBinning/DASTool/unbinned/'
        }
    }
    binning_depths_bins {
        path 'GenomeBinning/depths/bins'
    }
    binning_dereplication {
        path 'GenomeBinning/Dereplication'
    }
    binning_summary {
        path 'GenomeBinning'
    }
    binning_contig2bin {
        path 'GenomeBinning/contig_to_bin'
    }
    binning_bigmag {
        path 'GenomeBinning/BIgMAG'
    }
    binqc_busco {
        path { meta, _files -> "GenomeBinning/QC/BUSCO/${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}" }
        overwrite false
    }
    binqc_busco_db {
        path 'GenomeBinning/QC/BUSCO'
        // several tasks publish the same shared `busco_downloads/lineages/` directory
        overwrite false
    }
    binqc_checkm_db {
        path 'GenomeBinning/QC/CheckM/checkm_downloads'
        overwrite false
    }
    binqc_checkm {
        path 'GenomeBinning/QC/CheckM'
    }
    binqc_checkm2_db {
        path 'GenomeBinning/QC/CheckM2/checkm2_downloads'
        overwrite false
    }
    binqc_checkm2 {
        path 'GenomeBinning/QC/CheckM2'
    }
    binqc_gunc_db {
        path 'GenomeBinning/QC/GUNC'
    }
    binqc_gunc_raw {
        path { meta, _files -> "GenomeBinning/QC/GUNC/raw/${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}" }
    }
    binqc_gunc_checkm {
        path { meta, _files -> "GenomeBinning/QC/GUNC/checkmmerged/${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}" }
    }
    binqc_summaries {
        path 'GenomeBinning/QC'
    }
    binqc_pydamage_bins {
        path 'GenomeBinning/QC/summarise_pydamagebins'
    }
    binqc_tiara {
        path 'GenomeBinning/Tiara'
    }
    taxonomy_tiara {
        path 'Taxonomy/Tiara'
    }
    taxonomy_cat_db {
        path 'Taxonomy/CAT/db'
    }
    taxonomy_cat_summary {
        path 'Taxonomy/CAT'
    }
    taxonomy_cat_bins {
        path { meta, _files -> "Taxonomy/CAT/${meta.assembler}/${meta.binner}/${meta.id}/bins" }
    }
    taxonomy_cat_unbins {
        path { meta, _files -> "Taxonomy/CAT/${meta.assembler}/${meta.binner}/${meta.id}/unbins" }
    }
    taxonomy_gtdbtk {
        path { meta, _dir -> "Taxonomy/GTDB-Tk/${meta.assembler}/${meta.binner}/${meta.id}" }
    }
    taxonomy_gtdbtk_summary {
        path 'Taxonomy/GTDB-Tk'
    }
    annotation_prodigal {
        path { meta, _files -> "Annotation/Prodigal/${meta.assembler}/${meta.id}" }
    }
    annotation_prokka {
        path { meta, _files -> "Annotation/Prokka/${meta.assembler}" }
    }
    annotation_bakta {
        path { meta, _files -> "Annotation/Bakta/${meta.assembler}" }
    }
    annotation_mmseqs_db {
        path 'Annotation/mmseqs_db'
    }
    annotation_metaeuk {
        path { meta, _files -> "Annotation/MetaEuk/${meta.assembler}" }
    }
    multiqc {
        path 'multiqc'
    }
}
