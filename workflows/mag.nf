/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_mag_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BINNING_PREPARATION             } from '../subworkflows/local/binning_preparation/main'
include { BINNING                         } from '../subworkflows/local/binning/main'
include { BIN_QC                          } from '../subworkflows/local/bin_qc/main'
include { BINNING_REFINEMENT              } from '../subworkflows/local/binning_refinement/main'
include { VIRUS_IDENTIFICATION            } from '../subworkflows/local/virus_identification/main'
include { GTDBTK                          } from '../subworkflows/local/gtdbtk/main'
include { ANCIENT_DNA_ASSEMBLY_VALIDATION } from '../subworkflows/local/ancient_dna/main'
include { DOMAIN_CLASSIFICATION           } from '../subworkflows/local/domain_classification/main'
include { DEPTHS                          } from '../subworkflows/local/depths/main'
include { LONGREAD_PREPROCESSING          } from '../subworkflows/local/preprocessing_longread/main'
include { SHORTREAD_PREPROCESSING         } from '../subworkflows/local/preprocessing_shortread/main'
include { ASSEMBLY                        } from '../subworkflows/local/assembly/main'
include { CATPACK                         } from '../subworkflows/local/catpack/main'

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_ASSEMBLYINPUT  } from '../modules/nf-core/gunzip'
include { PRODIGAL                        } from '../modules/nf-core/prodigal/main'
include { PROKKA                          } from '../modules/nf-core/prokka/main'
include { MMSEQS_DATABASES                } from '../modules/nf-core/mmseqs/databases/main'
include { METAEUK_EASYPREDICT             } from '../modules/nf-core/metaeuk/easypredict/main'

//
// MODULE: Local to the pipeline
//
include { QUAST                           } from '../modules/local/quast_run/main'
include { QUAST_BINS                      } from '../modules/local/quast_bins/main'
include { QUAST_BINS_SUMMARY              } from '../modules/local/quast_bins_summary/main'
include { BIN_SUMMARY                     } from '../modules/local/bin_summary/main'

workflow MAG {
    take:
    ch_raw_short_reads  // channel: samplesheet read in from --input
    ch_raw_long_reads
    ch_input_assemblies

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ////////////////////////////////////////////////////
    /* --  Create channel for reference databases  -- */
    ////////////////////////////////////////////////////

    if (params.host_genome) {
        host_fasta = params.genomes[params.host_genome].fasta ?: false
        ch_host_fasta = Channel.value(file("${host_fasta}"))
        host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
        ch_host_bowtie2index = Channel.fromPath("${host_bowtie2index}", checkIfExists: true).first()
    }
    else if (params.host_fasta) {
        ch_host_fasta = Channel.fromPath("${params.host_fasta}", checkIfExists: true).first() ?: false

        if (params.host_fasta_bowtie2index) {
            ch_host_bowtie2index = Channel.fromPath("${params.host_fasta_bowtie2index}", checkIfExists: true).first()
        }
        else {
            ch_host_bowtie2index = Channel.empty()
        }
    }
    else {
        ch_host_fasta = Channel.empty()
        ch_host_bowtie2index = Channel.empty()
    }

    if (!params.keep_phix) {
        ch_phix_db_file = params.phix_reference ? Channel.value(file("${params.phix_reference}", checkIfExists: true)) : Channel.value(file("${projectDir}/assets/data/GCA_002596845.1_ASM259684v1_genomic.fna.gz", checkIfExists: true))
    }
    else {
        ch_phix_db_file = Channel.empty()
    }

    if (!params.keep_lambda) {
        ch_lambda_db = params.lambda_reference ? Channel.value(file("${params.lambda_reference}", checkIfExists: true)) : Channel.value(file("${projectDir}/assets/data/GCA_000840245.1_ViralProj14204_genomic.fna.gz", checkIfExists: true))
    }
    else {
        ch_lambda_db = Channel.value([])
    }

    if (params.genomad_db) {
        ch_genomad_db = file(params.genomad_db, checkIfExists: true)
    }
    else {
        ch_genomad_db = Channel.empty()
    }

    gtdb = params.skip_binqc || params.skip_gtdbtk ? false : params.gtdb_db

    if (gtdb) {
        gtdb = file("${gtdb}", checkIfExists: true)
    }
    else {
        gtdb = []
    }

    if (params.metaeuk_db && !params.skip_metaeuk) {
        ch_metaeuk_db = Channel.value(file("${params.metaeuk_db}", checkIfExists: true))
    }
    else {
        ch_metaeuk_db = Channel.empty()
    }

    // Get mmseqs db for MetaEuk if requested
    if (!params.skip_metaeuk && params.metaeuk_mmseqs_db) {
        MMSEQS_DATABASES(params.metaeuk_mmseqs_db)
        ch_versions = ch_versions.mix(MMSEQS_DATABASES.out.versions)
        ch_metaeuk_db = MMSEQS_DATABASES.out.database
    }

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    if (!params.assembly_input) {
        SHORTREAD_PREPROCESSING(
            ch_raw_short_reads,
            ch_host_fasta,
            ch_host_bowtie2index,
            ch_phix_db_file,
            params.skip_shortread_qc,
        )
        ch_versions = ch_versions.mix(SHORTREAD_PREPROCESSING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.multiqc_files.collect { it[1] }.ifEmpty([]))
        ch_short_reads = SHORTREAD_PREPROCESSING.out.short_reads
        ch_short_reads_assembly = SHORTREAD_PREPROCESSING.out.short_reads_assembly
    }
    else {
        ch_short_reads = ch_raw_short_reads.map { meta, reads ->
            def meta_new = meta - meta.subMap('run')
            [meta_new, reads]
        }
    }

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */

    LONGREAD_PREPROCESSING(
        ch_raw_long_reads,
        ch_short_reads,
        ch_lambda_db,
        ch_host_fasta,
        params.skip_longread_qc,
    )
    ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_PREPROCESSING.out.multiqc_files.collect { it[1] }.ifEmpty([]))
    ch_long_reads = LONGREAD_PREPROCESSING.out.long_reads

    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

    if (!params.assembly_input) {

        ASSEMBLY(
            ch_short_reads_assembly,
            ch_long_reads,
        )
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        ch_shortread_assemblies = ASSEMBLY.out.shortread_assemblies
        ch_longread_assemblies = ASSEMBLY.out.longread_assemblies

        ch_assemblies = ch_shortread_assemblies.mix(ch_longread_assemblies)
    }
    else {
        ch_assemblies_split = ch_input_assemblies.branch { _meta, assembly ->
            gzipped: assembly.getExtension() == "gz"
            ungzip: true
        }

        GUNZIP_ASSEMBLYINPUT(ch_assemblies_split.gzipped)
        ch_versions = ch_versions.mix(GUNZIP_ASSEMBLYINPUT.out.versions)

        ch_assemblies = Channel.empty()
        ch_assemblies = ch_assemblies.mix(ch_assemblies_split.ungzip, GUNZIP_ASSEMBLYINPUT.out.gunzip)
        ch_shortread_assemblies = ch_assemblies.filter { it[0].assembler.toUpperCase() in ['SPADES', 'SPADESHYBRID', 'MEGAHIT'] }
        ch_longread_assemblies = ch_assemblies.filter { it[0].assembler.toUpperCase() in ['FLYE', 'METAMDBG'] }
    }

    if (!params.skip_quast) {
        QUAST(ch_assemblies)
        ch_versions = ch_versions.mix(QUAST.out.versions)
    }

    /*
    ================================================================================
                                    Predict proteins
    ================================================================================
    */

    if (!params.skip_prodigal) {
        PRODIGAL(
            ch_assemblies,
            'gff',
        )
        ch_versions = ch_versions.mix(PRODIGAL.out.versions)
    }

    /*
    ================================================================================
                                    Virus identification
    ================================================================================
    */

    if (params.run_virus_identification) {
        VIRUS_IDENTIFICATION(ch_assemblies, ch_genomad_db)
        ch_versions = ch_versions.mix(VIRUS_IDENTIFICATION.out.versions)
    }

    /*
    ================================================================================
                                Binning preparation
    ================================================================================
    */

    if (!params.skip_binning || params.ancient_dna) {
        BINNING_PREPARATION(
            ch_shortread_assemblies,
            ch_short_reads,
            ch_longread_assemblies,
            ch_long_reads,
        )
        ch_versions = ch_versions.mix(BINNING_PREPARATION.out.versions)
    }

    /*
    ================================================================================
                                    Ancient DNA
    ================================================================================
    */

    if (params.ancient_dna) {
        ANCIENT_DNA_ASSEMBLY_VALIDATION(BINNING_PREPARATION.out.grouped_mappings)
        ch_versions = ch_versions.mix(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.versions)
    }

    /*
    ================================================================================
                                    Binning
    ================================================================================
    */

    if (!params.skip_binning) {

        // Make sure if running aDNA subworkflow to use the damage-corrected contigs for higher accuracy
        if (params.ancient_dna && !params.skip_ancient_damagecorrection) {
            BINNING(
                BINNING_PREPARATION.out.grouped_mappings.join(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled).map { it -> [it[0], it[4], it[2], it[3]] },
                params.bin_min_size,
                params.bin_max_size,
            )
        }
        else {
            BINNING(
                BINNING_PREPARATION.out.grouped_mappings,
                params.bin_min_size,
                params.bin_max_size,
            )
        }
        ch_versions = ch_versions.mix(BINNING.out.versions)

        if (params.bin_domain_classification) {

            // Make sure if running aDNA subworkflow to use the damage-corrected contigs for higher accuracy
            if (params.ancient_dna && !params.skip_ancient_damagecorrection) {
                ch_assemblies_for_domainclassification = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
            }
            else {
                ch_assemblies_for_domainclassification = ch_assemblies
            }

            DOMAIN_CLASSIFICATION(ch_assemblies_for_domainclassification, BINNING.out.bins, BINNING.out.unbinned)
            ch_versions = ch_versions.mix(DOMAIN_CLASSIFICATION.out.versions)

            ch_binning_results_bins = DOMAIN_CLASSIFICATION.out.classified_bins
            ch_binning_results_unbins = DOMAIN_CLASSIFICATION.out.classified_unbins
        }
        else {
            ch_binning_results_bins = BINNING.out.bins.map { meta, bins ->
                def meta_new = meta + [domain: 'unclassified']
                [meta_new, bins]
            }
            ch_binning_results_unbins = BINNING.out.unbinned.map { meta, bins ->
                def meta_new = meta + [domain: 'unclassified']
                [meta_new, bins]
            }
        }

        /*
        * DAS Tool: binning refinement
        */

        ch_binning_results_bins = ch_binning_results_bins.map { meta, bins ->
            def meta_new = meta + [refinement: 'unrefined']
            [meta_new, bins]
        }

        ch_binning_results_unbins = ch_binning_results_unbins.map { meta, bins ->
            def meta_new = meta + [refinement: 'unrefined_unbinned']
            [meta_new, bins]
        }

        // If any two of the binners are both skipped at once, do not run because DAS_Tool needs at least one
        if (params.refine_bins_dastool) {
            ch_prokarya_bins_dastool = ch_binning_results_bins.filter { meta, _bins ->
                meta.domain != "eukarya"
            }

            if (params.ancient_dna) {
                ch_contigs_for_binrefinement = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
            }
            else {
                ch_contigs_for_binrefinement = BINNING_PREPARATION.out.grouped_mappings.map { meta, contigs, _bam, _bai -> [meta, contigs] }
            }

            BINNING_REFINEMENT(ch_contigs_for_binrefinement, ch_prokarya_bins_dastool)
            ch_versions = ch_versions.mix(BINNING_REFINEMENT.out.versions)

            ch_refined_bins = BINNING_REFINEMENT.out.refined_bins
            ch_refined_unbins = BINNING_REFINEMENT.out.refined_unbins

            if (params.postbinning_input == 'raw_bins_only') {
                ch_input_for_postbinning_bins = ch_binning_results_bins
                ch_input_for_postbinning_unbins = ch_binning_results_unbins
            }
            else if (params.postbinning_input == 'refined_bins_only') {
                ch_input_for_postbinning_bins = ch_refined_bins
                ch_input_for_postbinning_unbins = ch_refined_unbins
            }
            else if (params.postbinning_input == 'both') {
                ch_input_for_postbinning_bins = ch_binning_results_bins.mix(ch_refined_bins)
                ch_input_for_postbinning_unbins = ch_binning_results_unbins.mix(ch_refined_unbins)
            }
        }
        else {
            ch_input_for_postbinning_bins = ch_binning_results_bins
            ch_input_for_postbinning_unbins = ch_binning_results_unbins
        }

        ch_input_for_postbinning = params.exclude_unbins_from_postbinning
            ? ch_input_for_postbinning_bins
            : ch_input_for_postbinning_bins.mix(ch_input_for_postbinning_unbins)

        // Combine short and long reads by meta.id and meta.group for DEPTHS, making sure that
        // read channel are not empty
        ch_reads_for_depths = ch_short_reads
            .map { meta, reads -> [[id: meta.id, group: meta.group], [short_reads: reads, long_reads: []]] }
            .mix(
                ch_long_reads.map { meta, reads -> [[id: meta.id, group: meta.group], [short_reads: [], long_reads: reads]] }
            )
            .groupTuple(by: 0)

        DEPTHS(ch_input_for_postbinning, BINNING.out.metabat2depths, ch_reads_for_depths)
        ch_versions = ch_versions.mix(DEPTHS.out.versions)

        ch_input_for_binsummary = DEPTHS.out.depths_summary

        /*
        * Bin QC subworkflows: for checking bin completeness with either BUSCO, CHECKM, CHECKM2, and/or GUNC
        */

        ch_bin_qc_summary = Channel.empty()
        if (!params.skip_binqc) {
            BIN_QC(ch_input_for_postbinning)
            ch_versions = ch_versions.mix(BIN_QC.out.versions)
            ch_bin_qc_summary = BIN_QC.out.qc_summaries
            ch_busco_summary = BIN_QC.out.busco_summary
            ch_checkm_summary = BIN_QC.out.checkm_summary
            ch_checkm2_summary = BIN_QC.out.checkm2_summary
        }

        ch_quast_bins_summary = Channel.empty()
        if (!params.skip_quast) {
            ch_input_for_quast_bins = ch_input_for_postbinning
                .groupTuple()
                .map { meta, bins ->
                    def new_bins = bins.flatten()
                    [meta, new_bins]
                }

            QUAST_BINS(ch_input_for_quast_bins)
            ch_versions = ch_versions.mix(QUAST_BINS.out.versions)
            ch_quast_bin_summary = QUAST_BINS.out.quast_bin_summaries.collectFile(keepHeader: true) { meta, summary ->
                ["${meta.id}.tsv", summary]
            }
            QUAST_BINS_SUMMARY(ch_quast_bin_summary.collect())
            ch_versions = ch_versions.mix(QUAST_BINS_SUMMARY.out.versions)

            ch_quast_bins_summary = QUAST_BINS_SUMMARY.out.summary
        }

        /*
         * CATPACK: bin / contig taxonomic classification
         */
        ch_catpack_summary = Channel.empty()
        if (params.cat_db || params.cat_db_generate) {
            CATPACK(
                ch_input_for_postbinning_bins,
                ch_input_for_postbinning_unbins,
            )
            ch_versions = ch_versions.mix(CATPACK.out.versions)

            ch_catpack_summary = CATPACK.out.summary
        }

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */

        if (!params.skip_gtdbtk) {
            ch_gtdbtk_summary = Channel.empty()
            if (gtdb) {

                ch_gtdb_bins = ch_input_for_postbinning.filter { meta, _bins ->
                    meta.domain != "eukarya"
                }

                GTDBTK(
                    ch_gtdb_bins,
                    ch_bin_qc_summary,
                    gtdb,
                )
                ch_versions = ch_versions.mix(GTDBTK.out.versions)
                ch_gtdbtk_summary = GTDBTK.out.summary
            }
        }
        else {
            ch_gtdbtk_summary = Channel.empty()
        }
        if ((!params.skip_binqc) || !params.skip_quast || !params.skip_gtdbtk) {
            BIN_SUMMARY(
                ch_input_for_binsummary,
                ch_quast_bins_summary.ifEmpty([]),
                ch_gtdbtk_summary.ifEmpty([]),
                ch_catpack_summary.ifEmpty([]),
                ch_busco_summary.ifEmpty([]),
                ch_checkm_summary.ifEmpty([]),
                ch_checkm2_summary.ifEmpty([]),
            )
            ch_versions = ch_versions.mix(BIN_SUMMARY.out.versions)
        }

        /*
         * Prokka: Genome annotation
         */

        if (!params.skip_prokka) {
            ch_bins_for_prokka = ch_input_for_postbinning
                .transpose()
                .map { meta, bin ->
                    def meta_new = meta + [id: bin.getBaseName()]
                    [meta_new, bin]
                }
                .filter { meta, _bin ->
                    meta.domain != "eukarya"
                }

            PROKKA(
                ch_bins_for_prokka,
                [],
                [],
            )
            ch_versions = ch_versions.mix(PROKKA.out.versions)
        }

        if (!params.skip_metaeuk && (params.metaeuk_db || params.metaeuk_mmseqs_db)) {
            ch_bins_for_metaeuk = ch_input_for_postbinning
                .transpose()
                .filter { meta, _bin ->
                    meta.domain in ["eukarya", "unclassified"]
                }
                .map { meta, bin ->
                    def meta_new = meta + [id: bin.getBaseName()]
                    [meta_new, bin]
                }

            METAEUK_EASYPREDICT(ch_bins_for_metaeuk, ch_metaeuk_db)
            ch_versions = ch_versions.mix(METAEUK_EASYPREDICT.out.versions)
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'mag_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.fromPath("${workflow.projectDir}/docs/images/mag_logo_mascot_light.png", checkIfExists: true)

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    if (!params.skip_quast) {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))

        if (!params.skip_binning) {
            ch_multiqc_files = ch_multiqc_files.mix(QUAST_BINS.out.dir.collect().ifEmpty([]))
        }
    }

    if (!params.skip_binning || params.ancient_dna) {
        ch_multiqc_files = ch_multiqc_files.mix(BINNING_PREPARATION.out.multiqc_files.collect().ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_prokka) {
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect { it[1] }.ifEmpty([]))
    }
    if (!params.skip_binning && !params.skip_binqc) {
        ch_multiqc_files = ch_multiqc_files.mix(BIN_QC.out.multiqc_files.collect().ifEmpty([]))
    }

    if (!params.skip_gtdbtk) {
        ch_multiqc_files = ch_multiqc_files.mix(GTDBTK.out.multiqc_files.collect().ifEmpty([]))
    }

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
