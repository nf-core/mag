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
    raw_short_reads  // channel: samplesheet read in from --input
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
    multiqc_report              = MAG.out.multiqc_report              // channel: /path/to/multiqc_report.html
    shortread_qc                = MAG.out.shortread_qc
    longread_qc                 = MAG.out.longread_qc
    assembly                    = MAG.out.assembly
    assembly_qc                 = MAG.out.assembly_qc
    virus_identification        = MAG.out.virus_identification
    ancient_dna                 = MAG.out.ancient_dna
    binning                     = MAG.out.binning
    bin_qc                      = MAG.out.bin_qc
    taxonomy                    = MAG.out.taxonomy
    annotation                  = MAG.out.annotation
    databases                   = MAG.out.databases
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
    shortread_qc                = NFCORE_MAG.out.shortread_qc
    longread_qc                 = NFCORE_MAG.out.longread_qc
    assembly                    = NFCORE_MAG.out.assembly
    assembly_qc                 = NFCORE_MAG.out.assembly_qc
    virus_identification        = NFCORE_MAG.out.virus_identification
    ancient_dna                 = NFCORE_MAG.out.ancient_dna
    binning                     = NFCORE_MAG.out.binning
    bin_qc                      = NFCORE_MAG.out.bin_qc
    taxonomy                    = NFCORE_MAG.out.taxonomy
    annotation                  = NFCORE_MAG.out.annotation
    databases                   = NFCORE_MAG.out.databases
    multiqc                     = NFCORE_MAG.out.multiqc
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

output {
    shortread_qc {
        path { meta, out ->
            out.fastqc >> 'QC_shortreads/fastqc/'
            out.clipping >> "QC_shortreads/${params.clip_tool}/${meta.id}/"
            out.clipping_reads >> (params.save_clipped_reads ? "QC_shortreads/${params.clip_tool}/${meta.id}/" : null)
            out.hostremoval >> 'QC_shortreads/remove_host/'
            out.hostremoval_reads >> (params.save_hostremoved_reads ? 'QC_shortreads/remove_host/' : null)
            out.phixremoval >> 'QC_shortreads/remove_phix/'
            out.phixremoval_reads >> (params.save_phixremoved_reads ? 'QC_shortreads/remove_phix/' : null)
            out.bbnorm_fastq >> (params.save_bbnorm_reads ? 'bbmap/bbnorm/' : null)
            out.bbnorm_log >> (params.save_bbnorm_reads ? 'bbmap/bbnorm/logs/' : null)
        }
    }
    longread_qc {
        path { meta, out ->
            out.adaptertrimming >> (params.save_porechop_reads ? 'QC_longreads/porechop/' : null)
            out.filtering >> "QC_longreads/${[filtlong: 'Filtlong', nanoq: 'Nanoq', chopper: 'Chopper'][params.longread_filtering_tool]}/"
            out.filtering_reads >> (params.save_filtered_longreads || (params.longread_filtering_tool == 'chopper' && params.save_lambdaremoved_reads) ? "QC_longreads/${[filtlong: 'Filtlong', nanoq: 'Nanoq', chopper: 'Chopper'][params.longread_filtering_tool]}/" : null)
            out.nanolyse >> 'QC_longreads/NanoLyse/'
            out.nanolyse_reads >> (params.save_lambdaremoved_reads ? 'QC_longreads/NanoLyse/' : null)
            out.nanoplot >> "QC_longreads/NanoPlot/${meta.id}/"
            out.hostremoval_index >> 'QC_longreads/minimap2/index/'
            out.hostremoval_bam >> (params.save_hostremoved_reads ? 'QC_longreads/minimap2/align/' : null)
            out.hostremoval_fastq >> (params.save_hostremoved_reads ? 'QC_longreads/samtools/fastq/' : null)
            out.hostremoval_stats >> 'QC_longreads/samtools/stats/'
        }
    }
    assembly {
        // pypolca vcf and report already sit in a <assembler>-<id>/ subdir of the task dir, which is kept
        path { meta, out ->
            out.assembly >> "Assembly/${meta.assembler}/"
            out.pypolca_polished >> "Assembly/${meta.assembler}-pypolca/${meta.assembler}-${meta.id}/${out.pypolca_polished?.name?.replace('_corrected.fasta', '_pypolca.fasta')}"
            out.pypolca_vcf >> "Assembly/${meta.assembler}-pypolca/"
            out.pypolca_report >> "Assembly/${meta.assembler}-pypolca/"
        }
    }
    assembly_qc {
        path { meta, out ->
            out.quast >> "Assembly/${meta.assembler}/QC/${meta.id}/"
            out.mapping >> "Assembly/${meta.assembler}/QC/${meta.id}/"
            out.mapping_bams >> (params.save_assembly_mapped_reads ? "Assembly/${meta.assembler}/QC/${meta.id}/" : null)
            out.mapping_full >> "Assembly/${meta.assembler}/QC/${meta.id}/full_assembly/"
            out.mapping_full_bams >> (params.save_assembly_mapped_reads ? "Assembly/${meta.assembler}/QC/${meta.id}/full_assembly/" : null)
            out.ale >> "Assembly/${meta.assembler}/QC/${meta.id}/ALE/"
            out.deepmased_features >> "Assembly/${meta.assembler}/QC/${meta.id}/DeepMAsED/features/"
            out.deepmased >> "Assembly/${meta.assembler}/QC/${meta.id}/DeepMAsED/"
        }
    }
    virus_identification {
        path { meta, out ->
            out.genomad >> "VirusIdentification/geNomad/${meta.id}/"
        }
    }
    ancient_dna {
        path { meta, out ->
            out.freebayes >> 'Ancient_DNA/variant_calling/freebayes/'
            out.filtered >> 'Ancient_DNA/variant_calling/filtered/'
            out.consensus >> 'Ancient_DNA/variant_calling/consensus/'
            out.pydamage_analyze >> "Ancient_DNA/pydamage/analyze/${meta.assembler}-${meta.id}/"
            out.pydamage_filter >> "Ancient_DNA/pydamage/filter/${meta.assembler}-${meta.id}/"
        }
    }
    binning {
        path { meta, out ->
            out.depths_contigs >> 'GenomeBinning/depths/contigs/'
            out.metabat2_bins >> 'GenomeBinning/MetaBAT2/bins/'
            out.metabat2_tooshort >> 'GenomeBinning/MetaBAT2/discarded/'
            out.metabat2_lowdepth >> 'GenomeBinning/MetaBAT2/discarded/'
            out.maxbin2_bins >> 'GenomeBinning/MaxBin2/bins/'
            out.maxbin2_tooshort >> 'GenomeBinning/MaxBin2/discarded/'
            out.maxbin2_summary >> 'GenomeBinning/MaxBin2/'
            out.maxbin2_abundance >> 'GenomeBinning/MaxBin2/'
            out.concoct_stats >> 'GenomeBinning/CONCOCT/stats/'
            out.concoct_bin >> "GenomeBinning/CONCOCT/bins/${out.concoct_bin?.name}"
            out.comebin_bins >> 'GenomeBinning/COMEBin/bins/'
            out.comebin_tsv >> 'GenomeBinning/COMEBin/stats/'
            out.comebin_log >> 'GenomeBinning/COMEBin/stats/'
            out.comebin_embeddings >> 'GenomeBinning/COMEBin/stats/'
            out.comebin_covembeddings >> 'GenomeBinning/COMEBin/stats/'
            out.metabinner_membership >> 'GenomeBinning/MetaBinner/stats/'
            out.metabinner_log >> 'GenomeBinning/MetaBinner/stats/'
            out.metabinner_bins >> 'GenomeBinning/MetaBinner/'
            out.metabinner_unbinned >> 'GenomeBinning/MetaBinner/unbinned/'
            out.semibin_bins >> 'GenomeBinning/SemiBin2/'
            out.semibin_log >> 'GenomeBinning/SemiBin2/log/'
            out.unbinned >> "GenomeBinning/${meta.binner}/unbinned/"
            out.unbinned_pooled >> "GenomeBinning/${meta.binner}/unbinned/discarded/"
            out.dastool_contig2bin >> (params.refine_bins_dastool_savecontig2bin ? 'GenomeBinning/DASTool/fastatocontig2bin/' : null)
            out.dastool >> 'GenomeBinning/DASTool/'
            out.dastool_bins >> 'GenomeBinning/DASTool/bins/'
            out.dastool_unbins >> 'GenomeBinning/DASTool/unbinned/'
            out.depths_bins >> 'GenomeBinning/depths/bins/'
            out.dereplication_tsv >> 'GenomeBinning/Dereplication/'
            out.dereplication_bins >> 'GenomeBinning/Dereplication/'
            out.summary >> 'GenomeBinning/'
            out.contig2bin >> 'GenomeBinning/contig_to_bin/'
            out.bigmag >> 'GenomeBinning/BIgMAG/'
        }
    }
    bin_qc {
        path { meta, out ->
            out.busco >> "GenomeBinning/QC/BUSCO/${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}/"
            out.checkm >> 'GenomeBinning/QC/CheckM/'
            out.checkm2 >> 'GenomeBinning/QC/CheckM2/'
            out.gunc_raw >> "GenomeBinning/QC/GUNC/raw/${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}/"
            out.gunc_checkm >> "GenomeBinning/QC/GUNC/checkmmerged/${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}/"
            out.summaries >> 'GenomeBinning/QC/'
            out.pydamage_bins >> 'GenomeBinning/QC/summarise_pydamagebins/'
            out.tiara >> 'GenomeBinning/Tiara/'
        }
    }
    taxonomy {
        path { meta, out ->
            out.tiara >> 'Taxonomy/Tiara/'
            out.cat_summary >> 'Taxonomy/CAT/'
            out.cat_bins >> "Taxonomy/CAT/${meta.assembler}/${meta.binner}/${meta.id}/bins/"
            out.cat_unbins >> "Taxonomy/CAT/${meta.assembler}/${meta.binner}/${meta.id}/unbins/"
            out.gtdbtk >> "Taxonomy/GTDB-Tk/${meta.assembler}/${meta.binner}/${meta.id}/"
            out.gtdbtk_summary >> 'Taxonomy/GTDB-Tk/'
        }
    }
    annotation {
        path { meta, out ->
            out.prodigal >> "Annotation/Prodigal/${meta.assembler}/${meta.id}/"
            out.prokka >> "Annotation/Prokka/${meta.assembler}/"
            out.bakta >> "Annotation/Bakta/${meta.assembler}/"
            out.metaeuk >> "Annotation/MetaEuk/${meta.assembler}/"
        }
    }
    databases {
        path { _meta, out ->
            out.busco_db >> (params.save_busco_db ? 'GenomeBinning/QC/BUSCO/' : null)
            out.checkm_db >> (params.save_checkm_data ? 'GenomeBinning/QC/CheckM/checkm_downloads/' : null)
            out.checkm2_db >> (params.save_checkm2_data ? 'GenomeBinning/QC/CheckM2/checkm2_downloads/' : null)
            out.gunc_db >> (params.gunc_save_db ? 'GenomeBinning/QC/GUNC/' : null)
            out.cat_db >> (params.save_cat_db ? 'Taxonomy/CAT/db/' : null)
            out.mmseqs_db >> (params.save_mmseqs_db ? 'Annotation/mmseqs_db/' : null)
        }
        // several BUSCO tasks publish the same shared `busco_downloads/lineages/` directory
        overwrite false
    }
    multiqc {
        path 'multiqc'
    }
}
