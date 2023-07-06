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
                if (row.long_reads) hybrid = true
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

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/nf-core-mag_logo_light.png")
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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
include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS        } from '../modules/local/pool_single_reads'
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
include { CAT_SUMMARY                                         } from "../modules/local/cat_summary"
include { BIN_SUMMARY                                         } from '../modules/local/bin_summary'
include { COMBINE_TSV as COMBINE_SUMMARY_TSV                  } from '../modules/local/combine_tsv'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { BINNING_PREPARATION } from '../subworkflows/local/binning_preparation'
include { BINNING             } from '../subworkflows/local/binning'
include { BINNING_REFINEMENT  } from '../subworkflows/local/binning_refinement'
include { BUSCO_QC            } from '../subworkflows/local/busco_qc'
include { CHECKM_QC           } from '../subworkflows/local/checkm_qc'
include { GUNC_QC             } from '../subworkflows/local/gunc_qc'
include { GTDBTK              } from '../subworkflows/local/gtdbtk'
include { ANCIENT_DNA_ASSEMBLY_VALIDATION } from '../subworkflows/local/ancient_dna'
include { DOMAIN_CLASSIFICATION } from '../subworkflows/local/domain_classification'
include { DEPTHS              } from '../subworkflows/local/depths'

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
include { FASTQC as FASTQC_TRIMMED               } from '../modules/nf-core/fastqc/main'
include { SEQTK_MERGEPE                          } from '../modules/nf-core/seqtk/mergepe/main'
include { BBMAP_BBNORM                           } from '../modules/nf-core/bbmap/bbnorm/main'
include { FASTP                                  } from '../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE    } from '../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE    } from '../modules/nf-core/adapterremoval/main'
include { CAT_FASTQ                              } from '../modules/nf-core/cat/fastq/main'
include { PRODIGAL                               } from '../modules/nf-core/prodigal/main'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if ( params.host_genome ) {
    host_fasta = params.genomes[params.host_genome].fasta ?: false
    ch_host_fasta = Channel
        .value(file( "${host_fasta}" ))
    host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
    ch_host_bowtie2index = Channel
        .value(file( "${host_bowtie2index}/*" ))
} else if ( params.host_fasta ) {
    ch_host_fasta = Channel
        .value(file( "${params.host_fasta}" ))
} else {
    ch_host_fasta = Channel.empty()
}

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

if(params.centrifuge_db){
    ch_centrifuge_db_file = Channel
        .value(file( "${params.centrifuge_db}" ))
} else {
    ch_centrifuge_db_file = Channel.empty()
}

if(params.kraken2_db){
    ch_kraken2_db_file = Channel
        .value(file( "${params.kraken2_db}" ))
} else {
    ch_kraken2_db_file = Channel.empty()
}

if(params.cat_db){
    ch_cat_db_file = Channel
        .value(file( "${params.cat_db}" ))
} else {
    ch_cat_db_file = Channel.empty()
}

if(!params.keep_phix) {
    ch_phix_db_file = Channel
        .value(file( "${params.phix_reference}" ))
}

if (!params.keep_lambda) {
    ch_nanolyse_db = Channel
        .value(file( "${params.lambda_reference}" ))
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

    ch_bowtie2_removal_host_multiqc = Channel.empty()
    if ( !params.assembly_input ) {
        if ( !params.skip_clipping ) {
            if ( params.clip_tool == 'fastp' ) {
                ch_clipmerge_out = FASTP (
                    ch_raw_short_reads,
                    [],
                    params.fastp_save_trimmed_fail,
                    []
                )
                ch_short_reads_prepped = FASTP.out.reads
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

            ch_short_reads_prepped = Channel.empty()
            ch_short_reads_prepped = ch_short_reads_prepped.mix(ADAPTERREMOVAL_SE.out.singles_truncated, ADAPTERREMOVAL_PE.out.paired_truncated)

                ch_versions = ch_versions.mix(ADAPTERREMOVAL_PE.out.versions.first(), ADAPTERREMOVAL_SE.out.versions.first())
            }
        } else {
            ch_short_reads_prepped = ch_raw_short_reads
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
                ch_short_reads_prepped,
                ch_host_bowtie2index
            )
            ch_short_reads_hostremoved = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
            ch_bowtie2_removal_host_multiqc = BOWTIE2_HOST_REMOVAL_ALIGN.out.log
            ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.versions.first())
        } else {
            ch_short_reads_hostremoved = ch_short_reads_prepped
        }

        if(!params.keep_phix) {
            BOWTIE2_PHIX_REMOVAL_BUILD (
                ch_phix_db_file
            )
            BOWTIE2_PHIX_REMOVAL_ALIGN (
                ch_short_reads_hostremoved,
                BOWTIE2_PHIX_REMOVAL_BUILD.out.index
            )
            ch_short_reads_phixremoved = BOWTIE2_PHIX_REMOVAL_ALIGN.out.reads
            ch_versions = ch_versions.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.versions.first())
        } else {
            ch_short_reads_phixremoved = ch_short_reads_hostremoved
        }

        if (!(params.keep_phix && params.skip_clipping && !(params.host_genome || params.host_fasta))) {
            FASTQC_TRIMMED (
                ch_short_reads_phixremoved
            )
            ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions)
        }

        // Run/Lane merging

        ch_short_reads_forcat = ch_short_reads_phixremoved
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                    meta_new.remove('run')
                [ meta_new, reads ]
            }
            .groupTuple()
            .branch {
                meta, reads ->
                    cat:       ( meta.single_end && reads.size() == 1 ) || ( !meta.single_end && reads.size() >= 2 )
                    skip_cat: true // Can skip merging if only single lanes
            }

        CAT_FASTQ ( ch_short_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ]} )

            ch_short_reads = Channel.empty()
            ch_short_reads = CAT_FASTQ.out.reads.mix( ch_short_reads_forcat.skip_cat ).map { meta, reads -> [ meta, reads.flatten() ]}
            ch_versions    = ch_versions.mix(CAT_FASTQ.out.versions.first())

        if ( params.bbnorm ) {
            if ( params.coassemble_group ) {
                // Interleave pairs, to be able to treat them as single ends when calling bbnorm. This prepares
                // for dropping the single_end parameter, but keeps assembly modules as they are, i.e. not
                // accepting a mix of single end and pairs.
                SEQTK_MERGEPE (
                    ch_short_reads.filter { ! it[0].single_end }
                )
                ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions.first())
                // Combine the interleaved pairs with any single end libraries. Set the meta.single_end to true (used by the bbnorm module).
                    ch_bbnorm = SEQTK_MERGEPE.out.reads
                        .mix(ch_short_reads.filter { it[0].single_end })
                        .map { [ [ id: sprintf("group%s", it[0].group), group: it[0].group, single_end: true ], it[1] ] }
                        .groupTuple()
            } else {
                ch_bbnorm = ch_short_reads
            }
            BBMAP_BBNORM ( ch_bbnorm )
            ch_versions = ch_versions.mix(BBMAP_BBNORM.out.versions)
            ch_short_reads_assembly = BBMAP_BBNORM.out.fastq
        } else {
            ch_short_reads_assembly = ch_short_reads
        }
    } else {
        ch_short_reads = ch_raw_short_reads
                        .map {
                            meta, reads ->
                                def meta_new = meta.clone()
                                meta_new.remove('run')
                            [ meta_new, reads ]
                        }
    }

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
                        .map {
                        meta, reads ->
                            def meta_new = meta.clone()
                            meta_new.remove('run')
                        [ meta_new, reads ]
                    }

    if ( !params.assembly_input ) {
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
        ch_short_reads_tmp = ch_short_reads
            .map { meta, sr -> [ meta.id, meta, sr ] }

        ch_short_and_long_reads = ch_long_reads
            .map { meta, lr -> [ meta.id, meta, lr ] }
            .join(ch_short_reads_tmp, by: 0)
            .map { id, meta_lr, lr, meta_sr, sr -> [ meta_lr, lr, sr[0], sr[1] ] }  // should not occur for single-end, since SPAdes (hybrid) does not support single-end

        FILTLONG (
            ch_short_and_long_reads
        )
        ch_long_reads = FILTLONG.out.reads
        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())

        NANOPLOT_FILTERED (
            ch_long_reads
        )
    }

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
        ch_tax_classifications = CENTRIFUGE.out.results_for_krona.mix(KRAKEN2.out.results_for_krona)
            . map { classifier, meta, report ->
                def meta_new = meta.clone()
                meta_new.classifier  = classifier
                [ meta_new, report ]
            }
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

    if ( !params.assembly_input ) {
        // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
        if (params.coassemble_group) {
            // short reads
            // group and set group as new id
            ch_short_reads_grouped = ch_short_reads_assembly
                .map { meta, reads -> [ meta.group, meta, reads ] }
                .groupTuple(by: 0)
                .map { group, metas, reads ->
                    def assemble_as_single = params.single_end || ( params.bbnorm && params.coassemble_group )
                    def meta         = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = assemble_as_single
                    if ( assemble_as_single ) [ meta, reads.collect { it }, [] ]
                    else [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                }
            // long reads
            // group and set group as new id
            ch_long_reads_grouped = ch_long_reads
                .map { meta, reads -> [ meta.group, meta, reads ] }
                .groupTuple(by: 0)
                .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    [ meta, reads.collect { it } ]
                }
        } else {
            ch_short_reads_grouped = ch_short_reads_assembly
                .filter { it[0].single_end }
                .map { meta, reads -> [ meta, [ reads ], [] ] }
                .mix (
                    ch_short_reads_assembly
                        .filter { ! it[0].single_end }
                        .map { meta, reads -> [ meta, [ reads[0] ], [ reads[1] ] ] }
                )
            ch_long_reads_grouped = ch_long_reads
        }

        ch_assemblies = Channel.empty()
        if (!params.skip_megahit){
            MEGAHIT ( ch_short_reads_grouped )
            ch_megahit_assemblies = MEGAHIT.out.assembly
                .map { meta, assembly ->
                    def meta_new = meta.clone()
                    meta_new.assembler  = "MEGAHIT"
                    [ meta_new, assembly ]
                }
            ch_assemblies = ch_assemblies.mix(ch_megahit_assemblies)
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        }

        // Co-assembly: pool reads for SPAdes
        if ( ! params.skip_spades || ! params.skip_spadeshybrid ){
            if ( params.coassemble_group ) {
                if ( params.bbnorm ) {
                    ch_short_reads_spades = ch_short_reads_grouped.map { [ it[0], it[1] ] }
                } else {
                    POOL_SHORT_SINGLE_READS (
                        ch_short_reads_grouped
                            .filter { it[0].single_end }
                    )
                    POOL_PAIRED_READS (
                        ch_short_reads_grouped
                            .filter { ! it[0].single_end }
                    )
                    ch_short_reads_spades = POOL_SHORT_SINGLE_READS.out.reads
                        .mix(POOL_PAIRED_READS.out.reads)
                }
            } else {
                ch_short_reads_spades = ch_short_reads_assembly
            }
            // long reads
            if (!params.single_end && !params.skip_spadeshybrid){
                POOL_LONG_READS ( ch_long_reads_grouped )
                ch_long_reads_spades = POOL_LONG_READS.out.reads
            } else {
                ch_long_reads_spades = Channel.empty()
            }
        } else {
            ch_short_reads_spades = Channel.empty()
            ch_long_reads_spades  = Channel.empty()
        }

        if (!params.single_end && !params.skip_spades){
            SPADES ( ch_short_reads_spades )
            ch_spades_assemblies = SPADES.out.assembly
                .map { meta, assembly ->
                    def meta_new = meta.clone()
                    meta_new.assembler  = "SPAdes"
                    [ meta_new, assembly ]
                }
            ch_assemblies = ch_assemblies.mix(ch_spades_assemblies)
            ch_versions = ch_versions.mix(SPADES.out.versions.first())
        }

        if (!params.single_end && !params.skip_spadeshybrid){
            ch_short_reads_spades_tmp = ch_short_reads_spades
                .map { meta, reads -> [ meta.id, meta, reads ] }

            ch_reads_spadeshybrid = ch_long_reads_spades
                .map { meta, reads -> [ meta.id, meta, reads ] }
                .combine(ch_short_reads_spades_tmp, by: 0)
                .map { id, meta_long, long_reads, meta_short, short_reads -> [ meta_short, long_reads, short_reads ] }

            SPADESHYBRID ( ch_reads_spadeshybrid )
            ch_spadeshybrid_assemblies = SPADESHYBRID.out.assembly
                .map { meta, assembly ->
                    def meta_new = meta.clone()
                    meta_new.assembler  = "SPAdesHybrid"
                    [ meta_new, assembly ]
                }
            ch_assemblies = ch_assemblies.mix(ch_spadeshybrid_assemblies)
            ch_versions = ch_versions.mix(SPADESHYBRID.out.versions.first())
        }
    } else {
        ch_assemblies = ch_input_assemblies
    }

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        QUAST ( ch_assemblies )
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

    ch_busco_summary            = Channel.empty()
    ch_checkm_summary           = Channel.empty()

    if ( !params.skip_binning || params.ancient_dna ) {
        BINNING_PREPARATION (
            ch_assemblies,
            ch_short_reads
        )
        ch_versions = ch_versions.mix(BINNING_PREPARATION.out.bowtie2_version.first())
    }

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

        // Make sure if running aDNA subworkflow to use the damage-corrected contigs for higher accuracy
        if (params.ancient_dna && !params.skip_ancient_damagecorrection) {
            BINNING (
                BINNING_PREPARATION.out.grouped_mappings
                    .join(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled)
                    .map { it -> [ it[0], it[4], it[2], it[3] ] }, // [meta, contigs_recalled, bam, bais]
                ch_short_reads
            )
        } else {
            BINNING (
                BINNING_PREPARATION.out.grouped_mappings,
                ch_short_reads
            )
        }
        ch_versions = ch_versions.mix(BINNING.out.versions)

        if ( params.bin_domain_classification ) {

            // Make sure if running aDNA subworkflow to use the damage-corrected contigs for higher accuracy
            if (params.ancient_dna && !params.skip_ancient_damagecorrection) {
                ch_assemblies_for_domainclassification = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
            } else {
                ch_assemblies_for_domainclassification = ch_assemblies
            }

            DOMAIN_CLASSIFICATION ( ch_assemblies_for_domainclassification, BINNING.out.bins, BINNING.out.unbinned )
            ch_binning_results_bins   = DOMAIN_CLASSIFICATION.out.classified_bins
            ch_binning_results_unbins = DOMAIN_CLASSIFICATION.out.classified_unbins
            ch_versions               = ch_versions.mix(DOMAIN_CLASSIFICATION.out.versions)


        } else {
            ch_binning_results_bins = BINNING.out.bins
                .map { meta, bins ->
                    meta_new = meta.clone()
                    meta_new.domain = 'unclassified'
                    [meta_new, bins]
                }
            ch_binning_results_unbins = BINNING.out.unbinned
                .map { meta, bins ->
                    meta_new = meta.clone()
                    meta_new.domain = 'unclassified'
                    [meta_new, bins]
                }
        }

        /*
        * DAS Tool: binning refinement
        */

        // If any two of the binners are both skipped at once, do not run because DAS_Tool needs at least one
        if ( params.refine_bins_dastool ) {
            ch_prokarya_bins_dastool = ch_binning_results_bins
                .filter { meta, bins ->
                    meta.domain != "eukarya"
                }

            ch_eukarya_bins_dastool = ch_binning_results_bins
                .filter { meta, bins ->
                    meta.domain == "eukarya"
                }

            if (params.ancient_dna) {
                ch_contigs_for_binrefinement = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
            } else {
                ch_contigs_for_binrefinement = BINNING_PREPARATION.out.grouped_mappings
                    .map{ meta, contigs, bam, bai -> [ meta, contigs ] }
            }

            BINNING_REFINEMENT ( ch_contigs_for_binrefinement, ch_prokarya_bins_dastool )
            ch_refined_bins = ch_eukarya_bins_dastool.mix(BINNING_REFINEMENT.out.refined_bins)
            ch_refined_unbins = BINNING_REFINEMENT.out.refined_unbins
            ch_versions = ch_versions.mix(BINNING_REFINEMENT.out.versions)

            if ( params.postbinning_input == 'raw_bins_only' ) {
                ch_input_for_postbinning_bins        = ch_binning_results_bins
                ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix(ch_binning_results_unbins)
            } else if ( params.postbinning_input == 'refined_bins_only' ) {
                ch_input_for_postbinning_bins        = ch_refined_bins
                ch_input_for_postbinning_bins_unbins = ch_refined_bins.mix(ch_refined_unbins)
            } else if ( params.postbinning_input == 'both' ) {
                ch_all_bins = ch_binning_results_bins.mix(ch_refined_bins)
                ch_input_for_postbinning_bins        = ch_all_bins
                ch_input_for_postbinning_bins_unbins = ch_all_bins.mix(ch_binning_results_unbins).mix(ch_refined_unbins)
            }
        } else {
            ch_input_for_postbinning_bins        = ch_binning_results_bins
            ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix(ch_binning_results_unbins)
        }

        DEPTHS ( ch_input_for_postbinning_bins_unbins, BINNING.out.metabat2depths, ch_short_reads )
        ch_input_for_binsummary = DEPTHS.out.depths_summary
        ch_versions = ch_versions.mix(DEPTHS.out.versions)

        /*
        * Bin QC subworkflows: for checking bin completeness with either BUSCO, CHECKM, and/or GUNC
        */

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

            ch_input_bins_for_checkm = ch_input_bins_for_qc
                .filter { meta, bins ->
                    meta.domain != "eukarya"
                }

            CHECKM_QC (
                ch_input_bins_for_checkm.groupTuple(),
                ch_checkm_db
            )
            ch_checkm_summary = CHECKM_QC.out.summary

            ch_versions = ch_versions.mix(CHECKM_QC.out.versions)

        }

        if ( params.run_gunc && params.binqc_tool == 'checkm' ) {
            GUNC_QC ( ch_input_bins_for_checkm, ch_gunc_db, CHECKM_QC.out.checkm_tsv )
            ch_versions = ch_versions.mix( GUNC_QC.out.versions )
        } else if ( params.run_gunc ) {
            ch_input_bins_for_gunc = ch_input_for_postbinning_bins_unbins
                .filter { meta, bins ->
                    meta.domain != "eukarya"
                }
            GUNC_QC ( ch_input_bins_for_qc, ch_gunc_db, [] )
            ch_versions = ch_versions.mix( GUNC_QC.out.versions )
        }

        ch_quast_bins_summary = Channel.empty()
        if (!params.skip_quast){
            ch_input_for_quast_bins = ch_input_for_postbinning_bins_unbins
                                        .groupTuple()
                                        .map {
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

            ch_gtdb_bins = ch_input_for_postbinning_bins_unbins
                .filter { meta, bins ->
                    meta.domain != "eukarya"
                }

            GTDBTK (
                ch_gtdb_bins,
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

        if (!params.skip_prokka){
            ch_bins_for_prokka = ch_input_for_postbinning_bins_unbins.transpose()
            .map { meta, bin ->
                def meta_new = meta.clone()
                meta_new.id  = bin.getBaseName()
                [ meta_new, bin ]
            }
            .filter { meta, bin ->
                meta.domain != "eukarya"
            }

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

    methods_description    = WorkflowMag.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))

    if (!params.assembly_input) {

        if ( !params.skip_clipping && params.clip_tool == 'adapterremoval' ) {
            ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_PE.out.settings.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_SE.out.settings.collect{it[1]}.ifEmpty([]))

        } else if ( !params.skip_clipping && params.clip_tool == 'fastp' )  {
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
        }

        if (!(params.keep_phix && params.skip_clipping && !(params.host_genome || params.host_fasta))) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))
        }

        if ( params.host_fasta || params.host_genome ) {
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
        }

        if(!params.keep_phix) {
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
        }

    }

    ch_multiqc_files = ch_multiqc_files.mix(CENTRIFUGE.out.kreport.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

    if (!params.skip_quast){
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))

        if ( !params.skip_binning ) {
            ch_multiqc_files = ch_multiqc_files.mix(QUAST_BINS.out.dir.collect().ifEmpty([]))
        }
    }

    if ( !params.skip_binning || params.ancient_dna ) {
        ch_multiqc_files = ch_multiqc_files.mix(BINNING_PREPARATION.out.bowtie2_assembly_multiqc.collect().ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_prokka){
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_binqc && params.binqc_tool == 'busco'){
        ch_multiqc_files = ch_multiqc_files.mix(BUSCO_QC.out.multiqc.collect().ifEmpty([]))
    }


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
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
