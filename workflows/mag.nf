/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check already if long reads are provided
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
def hybrid = false
if(hasExtension(params.input, "csv")){
    Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() == 5) {
                    if (row.long_reads) hybrid = true
                } else {
                    error("Input samplesheet contains row with ${row.size()} column(s). Expects 5.")
                }
            }
}

// Validate input parameters
WorkflowMag.initialise(params, log, hybrid)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.phix_reference, params.host_fasta, params.centrifuge_db, params.kraken2_db, params.cat_db, params.gtdb, params.lambda_reference, params.busco_reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
// Currently not used as using local version of MultiQC
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { QUAST                                               } from '../modules/local/quast'
include { QUAST_BINS                                          } from '../modules/local/quast_bins'
include { QUAST_BINS_SUMMARY                                  } from '../modules/local/quast_bins_summary'
include { CAT_DB                                              } from '../modules/local/cat_db'
include { CAT_DB_GENERATE                                     } from '../modules/local/cat_db_generate'
include { CAT                                                 } from '../modules/local/cat'
include { CAT_SUMMARY                                         } from "../modules/local/cat_summary"
include { BIN_SUMMARY                                         } from '../modules/local/bin_summary'
include { COMBINE_TSV as COMBINE_SUMMARY_TSV                  } from '../modules/local/combine_tsv'
include { MULTIQC                                             } from '../modules/local/multiqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                     } from '../subworkflows/local/input_check'
include { SHORT_READ_PREPROCESS           } from '../subworkflows/local/short_read_preprocess'
include { LONG_READ_PREPROCESS            } from '../subworkflows/local/long_read_preprocess'
include { SHORT_READ_TAXONOMY             } from '../subworkflows/local/short_read_taxonomy'
include { ASSEMBLY                        } from '../subworkflows/local/assembly'
include { BINNING_PREPARATION             } from '../subworkflows/local/binning_preparation'
include { BINNING                         } from '../subworkflows/local/binning'
include { BINNING_REFINEMENT              } from '../subworkflows/local/binning_refinement'
include { BUSCO_QC                        } from '../subworkflows/local/busco_qc'
include { CHECKM_QC                       } from '../subworkflows/local/checkm_qc'
include { GUNC_QC                         } from '../subworkflows/local/gunc_qc'
include { GTDBTK                          } from '../subworkflows/local/gtdbtk'
include { ANCIENT_DNA_ASSEMBLY_VALIDATION } from '../subworkflows/local/ancient_dna'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { ARIA2 as ARIA2_UNTAR                   } from '../modules/nf-core/aria2/main'
include { FASTQC as FASTQC_RAW                   } from '../modules/nf-core/fastqc/main'
include { PRODIGAL                               } from '../modules/nf-core/prodigal/main'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if(params.busco_reference){
    ch_busco_db_file = Channel
        .value(file( "${params.busco_reference}" ))
} else {
    ch_busco_db_file = Channel.empty()
}
if (params.busco_download_path) {
    ch_busco_download_folder = Channel
        .value(file( "${params.busco_download_path}" ))
} else {
    ch_busco_download_folder = Channel.empty()
}

if(params.checkm_db) {
    ch_checkm_db = file(params.checkm_db, checkIfExists: true)
}

if (params.gunc_db) {
    ch_gunc_db = file(params.gunc_db, checkIfExists: true)
} else {
    ch_gunc_db = Channel.empty()
}

if(params.cat_db){
    ch_cat_db_file = Channel
        .value(file( "${params.cat_db}" ))
} else {
    ch_cat_db_file = Channel.empty()
}

gtdb = params.skip_binqc ? false : params.gtdb
if (gtdb) {
    ch_gtdb = Channel
        .value(file( "${gtdb}" ))
} else {
    ch_gtdb = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report    = []
def busco_failed_bins = [:]

workflow MAG {

    ch_versions = Channel.empty()

    // Get checkM database if not supplied

    if ( !params.skip_binqc && params.binqc_tool == 'checkm' && !params.checkm_db ) {
        ARIA2_UNTAR (params.checkm_download_url)
        ch_checkm_db = ARIA2_UNTAR.out.downloaded_file
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_raw_short_reads  = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads   = INPUT_CHECK.out.raw_long_reads
    ch_input_assemblies = INPUT_CHECK.out.input_assemblies

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    FASTQC_RAW (
        ch_raw_short_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    ch_short_reads_preprocess = params.assembly_input ? Channel.empty() : ch_raw_short_reads

    SHORT_READ_PREPROCESS(ch_short_reads_preprocess)
    ch_versions = ch_versions.mix(SHORT_READ_PREPROCESS.out.versions)

    if ( !params.assembly_input ) {
        ch_short_reads            = SHORT_READ_PREPROCESS.out.short_reads
        ch_short_reads_assembly   = SHORT_READ_PREPROCESS.out.short_reads_assembly
    } else {
        ch_short_reads            = ch_raw_short_reads
        ch_short_reads_assembly   = Channel.empty()
    }

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */

    ch_long_reads_preprocess = params.assembly_input ? Channel.empty() : ch_raw_long_reads

    LONG_READ_PREPROCESS(ch_long_reads_preprocess, ch_short_reads)
    ch_versions = ch_versions.mix(LONG_READ_PREPROCESS.out.versions)

    if ( !params.assembly_input ) {
        ch_long_reads = LONG_READ_PREPROCESS.out.long_reads
    } else {
        ch_long_reads = Channel.empty()
    }
    /*
    ================================================================================
                                    Taxonomic information
    ================================================================================
    */

    SHORT_READ_TAXONOMY(ch_short_reads)
    ch_versions = ch_versions.mix(SHORT_READ_TAXONOMY.out.versions)

    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

    ch_assemblies = Channel.empty()
    ASSEMBLY(ch_short_reads_assembly, ch_long_reads)
    ch_assemblies = ch_assemblies.mix(ch_input_assemblies, ASSEMBLY.out.assemblies)
    ch_versions   = ch_versions.mix(ASSEMBLY.out.versions)

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        QUAST ( ch_assemblies )
        ch_quast_multiqc = QUAST.out.qc
        ch_versions      = ch_versions.mix(QUAST.out.versions.first())
    }

    /*
    ================================================================================
                                    Predict proteins
    ================================================================================
    */

    if (!params.skip_prodigal){
        PRODIGAL (
            ch_assemblies,
            'gff'
        )
        ch_versions = ch_versions.mix(PRODIGAL.out.versions.first())
    }

    /*
    ================================================================================
                                Binning preparation
    ================================================================================
    */


    ch_bowtie2_assembly_multiqc = Channel.empty()
    ch_busco_summary            = Channel.empty()
    ch_checkm_summary           = Channel.empty()
    ch_busco_multiqc            = Channel.empty()



    BINNING_PREPARATION (
        ch_assemblies,
        ch_short_reads
    )

    /*
    ================================================================================
                                    Ancient DNA
    ================================================================================
    */

    if (params.ancient_dna){
        ANCIENT_DNA_ASSEMBLY_VALIDATION(BINNING_PREPARATION.out.grouped_mappings)
        ch_versions = ch_versions.mix(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.versions.first())
    }

    /*
    ================================================================================
                                    Binning
    ================================================================================
    */

    if (!params.skip_binning){

        if (params.ancient_dna) {
            BINNING (
                BINNING_PREPARATION.out.grouped_mappings
                    .join(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled)
                    .map{ it -> [ it[0], it[4], it[2], it[3] ] }, // [meta, contigs_recalled, bam, bais]
                ch_short_reads
            )
        } else {
            BINNING (
                BINNING_PREPARATION.out.grouped_mappings,
                ch_short_reads
            )
        }

        ch_bowtie2_assembly_multiqc = BINNING_PREPARATION.out.bowtie2_assembly_multiqc
        ch_versions = ch_versions.mix(BINNING_PREPARATION.out.bowtie2_version.first())
        ch_versions = ch_versions.mix(BINNING.out.versions)

        /*
        * DAS Tool: binning refinement
        */

        // If any two of the binners are both skipped at once, do not run because DAS_Tool needs at least one
        if ( params.refine_bins_dastool ) {

            BINNING_REFINEMENT ( BINNING_PREPARATION.out.grouped_mappings, BINNING.out.bins, BINNING.out.metabat2depths, ch_short_reads )
            ch_versions = ch_versions.mix(BINNING_REFINEMENT.out.versions)

            if ( params.postbinning_input == 'raw_bins_only' ) {
                ch_input_for_postbinning_bins        = BINNING.out.bins
                ch_input_for_postbinning_bins_unbins = BINNING.out.bins.mix(BINNING.out.unbinned)
                ch_input_for_binsummary              = BINNING.out.depths_summary
            } else if ( params.postbinning_input == 'refined_bins_only' ) {
                ch_input_for_postbinning_bins        = BINNING_REFINEMENT.out.refined_bins
                ch_input_for_postbinning_bins_unbins = BINNING_REFINEMENT.out.refined_bins.mix(BINNING_REFINEMENT.out.refined_unbins)
                ch_input_for_binsummary              = BINNING_REFINEMENT.out.refined_depths_summary
            } else if (params.postbinning_input == 'both') {
                ch_input_for_postbinning_bins        = BINNING.out.bins.mix(BINNING_REFINEMENT.out.refined_bins)
                ch_input_for_postbinning_bins_unbins = BINNING.out.bins.mix(BINNING.out.unbinned,BINNING_REFINEMENT.out.refined_bins,BINNING_REFINEMENT.out.refined_unbins)
                ch_combinedepthtsvs_for_binsummary   = BINNING.out.depths_summary.mix(BINNING_REFINEMENT.out.refined_depths_summary)
                ch_input_for_binsummary              = COMBINE_SUMMARY_TSV ( ch_combinedepthtsvs_for_binsummary.collect() ).combined
            }

        } else {
                ch_input_for_postbinning_bins        = BINNING.out.bins
                ch_input_for_postbinning_bins_unbins = BINNING.out.bins.mix(BINNING.out.unbinned)
                ch_input_for_binsummary              = BINNING.out.depths_summary
        }

        /*
        * Bin QC subworkflows: for checking bin completeness with either BUSCO, CHECKM, and/or GUNC
        */

        // Results in: [ [meta], path_to_bin.fa ]
        ch_input_bins_for_qc = ch_input_for_postbinning_bins_unbins.transpose()

        if (!params.skip_binqc && params.binqc_tool == 'busco'){
            /*
            * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
            */
            BUSCO_QC (
                ch_busco_db_file,
                ch_busco_download_folder,
                ch_input_bins_for_qc
            )
            ch_busco_summary = BUSCO_QC.out.summary
            ch_busco_multiqc = BUSCO_QC.out.multiqc
            ch_versions = ch_versions.mix(BUSCO_QC.out.versions.first())
            // process information if BUSCO analysis failed for individual bins due to no matching genes
            BUSCO_QC.out
                .failed_bin
                .splitCsv(sep: '\t')
                .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
        }

        if (!params.skip_binqc && params.binqc_tool == 'checkm'){
            /*
            * CheckM subworkflow: Quantitative measures for the assessment of genome assembly
            */
            CHECKM_QC (
                ch_input_bins_for_qc.groupTuple(),
                ch_checkm_db
            )
            ch_checkm_summary = CHECKM_QC.out.summary

            // TODO custom output parsing? Add to MultiQC?
            ch_versions = ch_versions.mix(CHECKM_QC.out.versions)

        }

        if ( params.run_gunc && params.binqc_tool == 'checkm' ) {
            GUNC_QC ( ch_input_bins_for_qc, ch_gunc_db, CHECKM_QC.out.checkm_tsv )
            ch_versions = ch_versions.mix( GUNC_QC.out.versions )
        } else if ( params.run_gunc ) {
            GUNC_QC ( ch_input_bins_for_qc, ch_gunc_db, [] )
            ch_versions = ch_versions.mix( GUNC_QC.out.versions )
        }

        ch_quast_bins_summary = Channel.empty()
        if (!params.skip_quast){
            ch_input_for_quast_bins = ch_input_for_postbinning_bins_unbins
                                        .groupTuple()
                                        .map{
                                            meta, reads ->
                                                def new_reads = reads.flatten()
                                                [meta, new_reads]
                                            }
            QUAST_BINS ( ch_input_for_quast_bins )
            ch_versions = ch_versions.mix(QUAST_BINS.out.versions.first())
            QUAST_BINS_SUMMARY ( QUAST_BINS.out.quast_bin_summaries.collect() )
            ch_quast_bins_summary = QUAST_BINS_SUMMARY.out.summary
        }

        /*
         * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
         */
        ch_cat_db = Channel.empty()
        if (params.cat_db){
            CAT_DB ( ch_cat_db_file )
            ch_cat_db = CAT_DB.out.db
        } else if (params.cat_db_generate){
            CAT_DB_GENERATE ()
            ch_cat_db = CAT_DB_GENERATE.out.db
        }
        CAT (
            ch_input_for_postbinning_bins_unbins,
            ch_cat_db
        )
        CAT_SUMMARY(
            CAT.out.tax_classification.collect()
        )
        ch_versions = ch_versions.mix(CAT.out.versions.first())
        ch_versions = ch_versions.mix(CAT_SUMMARY.out.versions)

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */
        ch_gtdbtk_summary = Channel.empty()
        if ( gtdb ){
            GTDBTK (
                ch_input_for_postbinning_bins_unbins,
                ch_busco_summary,
                ch_checkm_summary,
                ch_gtdb
            )
            ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
            ch_gtdbtk_summary = GTDBTK.out.summary
        }

        if ( ( !params.skip_binqc ) || !params.skip_quast || gtdb){
            BIN_SUMMARY (
                ch_input_for_binsummary,
                ch_busco_summary.ifEmpty([]),
                ch_checkm_summary.ifEmpty([]),
                ch_quast_bins_summary.ifEmpty([]),
                ch_gtdbtk_summary.ifEmpty([])
            )
        }

        /*
         * Prokka: Genome annotation
         */
        ch_bins_for_prokka = ch_input_for_postbinning_bins_unbins.transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.id  = bin.getBaseName()
                [ meta_new, bin ]
            }

        if (!params.skip_prokka){
            PROKKA (
                ch_bins_for_prokka,
                [],
                []
            )
            ch_versions = ch_versions.mix(PROKKA.out.versions.first())
        }
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMag.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    // Currently not used due to local MultiQC module
    //methods_description    = WorkflowMag.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    //ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    /* //This is the template input with the nf-core module
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bowtie2_removal_host_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_quast_multiqc.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bowtie2_assembly_multiqc.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_busco_multiqc.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    */

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
        SHORT_READ_PREPROCESS.out.multiqc_fastqc_trimmed.collect().ifEmpty([]),
        SHORT_READ_PREPROCESS.out.multiqc_bowtie2_removal_host.collect{it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect().ifEmpty([]),
        ch_bowtie2_assembly_multiqc.collect().ifEmpty([]),
        ch_busco_multiqc.collect().ifEmpty([]),
        SHORT_READ_PREPROCESS.out.multiqc_readprep.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, busco_failed_bins)
    }
    NfcoreTemplate.summary(workflow, params, log, busco_failed_bins)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
