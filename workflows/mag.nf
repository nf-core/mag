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
                    log.error "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    System.exit(1)
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

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'
include { PORECHOP                                            } from '../modules/local/porechop'
include { NANOLYSE                                            } from '../modules/local/nanolyse'
include { FILTLONG                                            } from '../modules/local/filtlong'
include { NANOPLOT as NANOPLOT_RAW                            } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_FILTERED                       } from '../modules/local/nanoplot'
include { CENTRIFUGE_DB_PREPARATION                           } from '../modules/local/centrifuge_db_preparation'
include { CENTRIFUGE                                          } from '../modules/local/centrifuge'
include { KRAKEN2_DB_PREPARATION                              } from '../modules/local/kraken2_db_preparation'
include { KRAKEN2                                             } from '../modules/local/kraken2'
include { KRONA_DB                                            } from '../modules/local/krona_db'
include { KRONA                                               } from '../modules/local/krona'
include { POOL_SINGLE_READS                                   } from '../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                   } from '../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                } from '../modules/local/pool_single_reads'
include { MEGAHIT                                             } from '../modules/local/megahit'
include { SPADES                                              } from '../modules/local/spades'
include { SPADESHYBRID                                        } from '../modules/local/spadeshybrid'
include { QUAST                                               } from '../modules/local/quast'
include { QUAST_BINS                                          } from '../modules/local/quast_bins'
include { QUAST_BINS_SUMMARY                                  } from '../modules/local/quast_bins_summary'
include { CAT_DB                                              } from '../modules/local/cat_db'
include { CAT_DB_GENERATE                                     } from '../modules/local/cat_db_generate'
include { CAT                                                 } from '../modules/local/cat'
include { BIN_SUMMARY                                         } from '../modules/local/bin_summary'
include { MULTIQC                                             } from '../modules/local/multiqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { BINNING_PREPARATION } from '../subworkflows/local/binning_preparation'
include { BINNING             } from '../subworkflows/local/binning'
include { BUSCO_QC            } from '../subworkflows/local/busco_qc'
include { GTDBTK              } from '../subworkflows/local/gtdbtk'
include { ANCIENT_DNA_ASSEMLY_VALIDATION } from '../subworkflows/local/ancient_dna'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW                   } from '../modules/nf-core/modules/fastqc/main'
include { FASTQC as FASTQC_TRIMMED               } from '../modules/nf-core/modules/fastqc/main'
include { FASTP                                  } from '../modules/nf-core/modules/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE    } from '../modules/nf-core/modules/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE    } from '../modules/nf-core/modules/adapterremoval/main'
include { PRODIGAL                               } from '../modules/nf-core/modules/prodigal/main'
include { PROKKA                                 } from '../modules/nf-core/modules/prokka/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if ( params.host_genome ) {
    host_fasta = params.genomes[params.host_genome].fasta ?: false
    Channel
        .value(file( "${host_fasta}" ))
        .set { ch_host_fasta }

    host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
    Channel
        .value(file( "${host_bowtie2index}/*" ))
        .set { ch_host_bowtie2index }
} else if ( params.host_fasta ) {
    Channel
        .value(file( "${params.host_fasta}" ))
        .set { ch_host_fasta }
} else {
    ch_host_fasta = Channel.empty()
}

if(params.busco_reference){
    Channel
        .value(file( "${params.busco_reference}" ))
        .set { ch_busco_db_file }
} else {
    ch_busco_db_file = Channel.empty()
}
if (params.busco_download_path) {
    Channel
        .value(file( "${params.busco_download_path}" ))
        .set { ch_busco_download_folder }
} else {
    ch_busco_download_folder = Channel.empty()
}

if(params.centrifuge_db){
    Channel
        .value(file( "${params.centrifuge_db}" ))
        .set { ch_centrifuge_db_file }
} else {
    ch_centrifuge_db_file = Channel.empty()
}

if(params.kraken2_db){
    Channel
        .value(file( "${params.kraken2_db}" ))
        .set { ch_kraken2_db_file }
} else {
    ch_kraken2_db_file = Channel.empty()
}

if(params.cat_db){
    Channel
        .value(file( "${params.cat_db}" ))
        .set { ch_cat_db_file }
} else {
    ch_cat_db_file = Channel.empty()
}

if(!params.keep_phix) {
    Channel
        .value(file( "${params.phix_reference}" ))
        .set { ch_phix_db_file }
}

if (!params.keep_lambda) {
    Channel
        .value(file( "${params.lambda_reference}" ))
        .set { ch_nanolyse_db }
}

gtdb = params.skip_busco ? false : params.gtdb
if (gtdb) {
    Channel
        .value(file( "${gtdb}" ))
        .set { ch_gtdb }
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

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads  = INPUT_CHECK.out.raw_long_reads

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    FASTQC_RAW (
        ch_raw_short_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    if ( params.clip_tool == 'fastp' ) {
        ch_clipmerge_out = FASTP (
            ch_raw_short_reads,
            params.fastp_save_trimmed_fail,
            []
        )
        ch_short_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

    } else if ( params.clip_tool == 'adapterremoval' ) {

        // due to strange output file scheme in AR2, have to manually separate
        // SE/PE to allow correct pulling of reads after.
        ch_adapterremoval_in = ch_raw_short_reads
            .branch {
                    single: it[0]['single_end']
                    paired: !it[0]['single_end']
                }

        ADAPTERREMOVAL_PE ( ch_adapterremoval_in.paired, [] )
        ADAPTERREMOVAL_SE ( ch_adapterremoval_in.single, [] )

        // pair1 and 2 come for PE data from separate output channels, so bring
        // this back together again here
        ch_adapterremoval_pe_out = Channel.empty()
        ch_adapterremoval_pe_out = ADAPTERREMOVAL_PE.out.pair1_truncated
            .join(ADAPTERREMOVAL_PE.out.pair2_truncated)
            .map {
                [ it[0], [it[1], it[2]] ]
            }

        ch_short_reads = Channel.empty()
        ch_short_reads = ch_short_reads.mix(ADAPTERREMOVAL_SE.out.singles_truncated, ch_adapterremoval_pe_out)

        ch_versions = ch_versions.mix(ADAPTERREMOVAL_PE.out.versions, ADAPTERREMOVAL_SE.out.versions)

    }

    if (params.host_fasta){
        BOWTIE2_HOST_REMOVAL_BUILD (
            ch_host_fasta
        )
        ch_host_bowtie2index = BOWTIE2_HOST_REMOVAL_BUILD.out.index
    }
    ch_bowtie2_removal_host_multiqc = Channel.empty()
    if (params.host_fasta || params.host_genome){
        BOWTIE2_HOST_REMOVAL_ALIGN (
            ch_short_reads,
            ch_host_bowtie2index
        )
        ch_short_reads = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
        ch_bowtie2_removal_host_multiqc = BOWTIE2_HOST_REMOVAL_ALIGN.out.log
        ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.versions.first())
    }

    if(!params.keep_phix) {
        BOWTIE2_PHIX_REMOVAL_BUILD (
            ch_phix_db_file
        )
        BOWTIE2_PHIX_REMOVAL_ALIGN (
            ch_short_reads,
            BOWTIE2_PHIX_REMOVAL_BUILD.out.index
        )
        ch_short_reads = BOWTIE2_PHIX_REMOVAL_ALIGN.out.reads
    }

    FASTQC_TRIMMED (
        ch_short_reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions)

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */
    NANOPLOT_RAW (
        ch_raw_long_reads
    )
    ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first())

    ch_long_reads = ch_raw_long_reads
    if (!params.skip_adapter_trimming) {
        PORECHOP (
            ch_raw_long_reads
        )
        ch_long_reads = PORECHOP.out.reads
        ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
    }

    if (!params.keep_lambda) {
        NANOLYSE (
            ch_long_reads,
            ch_nanolyse_db
        )
        ch_long_reads = NANOLYSE.out.reads
        ch_versions = ch_versions.mix(NANOLYSE.out.versions.first())
    }

    // join long and short reads by sample name
    ch_short_reads
        .map { meta, sr -> [ meta.id, meta, sr ] }
        .set {ch_short_reads_tmp}

    ch_long_reads
        .map { meta, lr -> [ meta.id, meta, lr ] }
        .join(ch_short_reads_tmp, by: 0)
        .map { id, meta_lr, lr, meta_sr, sr -> [ meta_lr, lr, sr[0], sr[1] ] }  // should not occur for single-end, since SPAdes (hybrid) does not support single-end
        .set{ ch_short_and_long_reads }

    FILTLONG (
        ch_short_and_long_reads
    )
    ch_long_reads = FILTLONG.out.reads
    ch_versions = ch_versions.mix(FILTLONG.out.versions.first())

    NANOPLOT_FILTERED (
        ch_long_reads
    )

    /*
    ================================================================================
                                    Taxonomic information
    ================================================================================
    */
    CENTRIFUGE_DB_PREPARATION ( ch_centrifuge_db_file )
    CENTRIFUGE (
        ch_short_reads,
        CENTRIFUGE_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(CENTRIFUGE.out.versions.first())

    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db_file
    )
    KRAKEN2 (
        ch_short_reads,
        KRAKEN2_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

    if (( params.centrifuge_db || params.kraken2_db ) && !params.skip_krona){
        KRONA_DB ()
        CENTRIFUGE.out.results_for_krona.mix(KRAKEN2.out.results_for_krona)
            . map { classifier, meta, report ->
                def meta_new = meta.clone()
                meta_new.classifier  = classifier
                [ meta_new, report ]
            }
            .set { ch_tax_classifications }
        KRONA (
            ch_tax_classifications,
            KRONA_DB.out.db.collect()
        )
        ch_versions = ch_versions.mix(KRONA.out.versions.first())
    }

    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

    // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = params.single_end
                    if (!params.single_end) [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                    else [ meta, reads.collect { it }, [] ]
            }
            .set { ch_short_reads_grouped }
        // long reads
        // group and set group as new id
        ch_long_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def meta = [:]
                meta.id          = "group-$group"
                meta.group       = group
                [ meta, reads.collect { it } ]
            }
            .set { ch_long_reads_grouped }
    } else {
        ch_short_reads
            .map { meta, reads ->
                    if (!params.single_end){ [ meta, [reads[0]], [reads[1]] ] }
                    else [ meta, [reads], [] ] }
            .set { ch_short_reads_grouped }
    }

    ch_assemblies = Channel.empty()
    if (!params.skip_megahit){
        MEGAHIT ( ch_short_reads_grouped )
        MEGAHIT.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "MEGAHIT"
                [ meta_new, assembly ]
            }
            .set { ch_megahit_assemblies }
        ch_assemblies = ch_assemblies.mix(ch_megahit_assemblies)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
    }

    // Co-assembly: pool reads for SPAdes
    if (params.coassemble_group) {
        // short reads
        if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)){
            if (params.single_end){
                POOL_SINGLE_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_SINGLE_READS.out.reads
            } else {
                POOL_PAIRED_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_PAIRED_READS.out.reads
            }
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid){
            POOL_LONG_READS ( ch_long_reads_grouped )
            ch_long_reads_spades = POOL_LONG_READS.out.reads
        }
    } else {
        ch_short_reads_spades = ch_short_reads
        ch_long_reads
            .map { meta, reads -> [ meta, [reads] ] }
            .set { ch_long_reads_spades }
    }

    if (!params.single_end && !params.skip_spades){
        SPADES ( ch_short_reads_spades )
        SPADES.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, assembly ]
            }
            .set { ch_spades_assemblies }
        ch_assemblies = ch_assemblies.mix(ch_spades_assemblies)
        ch_versions = ch_versions.mix(SPADES.out.versions.first())
    }

    if (!params.single_end && !params.skip_spadeshybrid){
        ch_short_reads_spades
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .set {ch_short_reads_spades_tmp}
        ch_long_reads_spades
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .combine(ch_short_reads_spades_tmp, by: 0)
            .map { id, meta_long, long_reads, meta_short, short_reads -> [ meta_short, long_reads, short_reads ] }
            .set { ch_reads_spadeshybrid }
        SPADESHYBRID ( ch_reads_spadeshybrid )
        SPADESHYBRID.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdesHybrid"
                [ meta_new, assembly ]
            }
            .set { ch_spadeshybrid_assemblies }
        ch_assemblies = ch_assemblies.mix(ch_spadeshybrid_assemblies)
        ch_versions = ch_versions.mix(SPADESHYBRID.out.versions.first())
    }

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        QUAST ( ch_assemblies )
        ch_quast_multiqc = QUAST.out.qc
        ch_versions = ch_versions.mix(QUAST.out.versions.first())
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
        ANCIENT_DNA_ASSEMLY_VALIDATION(BINNING_PREPARATION.out.grouped_mappings)
        ch_versions = ch_versions.mix(ANCIENT_DNA_ASSEMLY_VALIDATION.out.versions.first())
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
                    .join(ANCIENT_DNA_ASSEMLY_VALIDATION.out.contigs_recalled)
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

        if (!params.skip_busco){
            /*
            * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
            */
            ch_input_bins_busco = BINNING.out.bins.mix( BINNING.out.unbinned ).transpose()
            BUSCO_QC (
                ch_busco_db_file,
                ch_busco_download_folder,
                ch_input_bins_busco
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

        ch_quast_bins_summary = Channel.empty()
        if (!params.skip_quast){
            ch_input_for_quast_bins = BINNING.out.bins
                                        .mix( BINNING.out.unbinned )
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
            BINNING.out.bins,
            ch_cat_db
        )
        ch_versions = ch_versions.mix(CAT.out.versions.first())

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */
        ch_gtdbtk_summary = Channel.empty()
        if ( gtdb ){
            GTDBTK (
                BINNING.out.bins,
                ch_busco_summary,
                ch_gtdb
            )
            ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
            ch_gtdbtk_summary = GTDBTK.out.summary
        }

        if (!params.skip_busco || !params.skip_quast || gtdb){
            BIN_SUMMARY (
                BINNING.out.depths_summary,
                ch_busco_summary.ifEmpty([]),
                ch_quast_bins_summary.ifEmpty([]),
                ch_gtdbtk_summary.ifEmpty([])
            )
        }

        /*
         * Prokka: Genome annotation
         */
        BINNING.out.bins.mix(BINNING.out.unbinned).transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.id  = bin.getBaseName()
                [ meta_new, bin ]
            }
            .set { ch_bins_for_prokka }

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
        ch_multiqc_files.collect()
    )
    */

    //This is the local module because (1) fastq file name clashes [probably solveable] and (2) host removal bowtie reporting
    ch_multiqc_additional = Channel.empty()

    if ( params.clip_tool == "fastp") {
        ch_multiqc_additional = ch_multiqc_additional.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    } else if ( params.clip_tool == "adapterremoval" ) {
        ch_multiqc_additional = ch_multiqc_additional.mix(ADAPTERREMOVAL_PE.out.log.collect{it[1]}.ifEmpty([]), ADAPTERREMOVAL_SE.out.log.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
        FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]),
        ch_bowtie2_removal_host_multiqc.collect{it[1]}.ifEmpty([]),
        ch_quast_multiqc.collect().ifEmpty([]),
        ch_bowtie2_assembly_multiqc.collect().ifEmpty([]),
        ch_busco_multiqc.collect().ifEmpty([]),
        ch_multiqc_additional.collect().ifEmpty([]),
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
