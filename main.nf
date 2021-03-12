#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mag
========================================================================================

nf-core/mag Analysis Pipeline. Started 2018-05-22.
#### Homepage / Documentation
https://github.com/nf-core/mag
#### Authors
Hadrien Gourlé HadrienG <hadrien.gourle@slu.se> - hadriengourle.com>
Daniel Straub <d4straub@gmail.com>
Sabrina Krakau <sabrinakrakau@gmail.com>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/mag --input 'samplesheet.csv' -profile docker"
    // nextflow run nf-core/mag --input '*_R{1,2}.fastq.gz' -profile docker
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* -- IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS-- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Local: Modules
include { GET_SOFTWARE_VERSIONS                               } from './modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']]          )
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from './modules/local/bowtie2_removal_build'       addParams( options: [:]                                   )
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from './modules/local/bowtie2_removal_align'       addParams( options: modules['bowtie2_host_removal_align'] )
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD } from './modules/local/bowtie2_removal_build'       addParams( options: [:]                                   )
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN } from './modules/local/bowtie2_removal_align'       addParams( options: modules['bowtie2_phix_removal_align'] )
include { PORECHOP                                            } from './modules/local/porechop'                    addParams( options: [:]                                   )
include { NANOLYSE                                            } from './modules/local/nanolyse'                    addParams( options: modules['nanolyse']                   )
include { FILTLONG                                            } from './modules/local/filtlong'                    addParams( options: [:]                                   )
include { NANOPLOT as NANOPLOT_RAW                            } from './modules/local/nanoplot'                    addParams( options: modules['nanoplot_raw']               )
include { NANOPLOT as NANOPLOT_FILTERED                       } from './modules/local/nanoplot'                    addParams( options: modules['nanoplot_filtered']          )
include { CENTRIFUGE_DB_PREPARATION                           } from './modules/local/centrifuge_db_preparation'   addParams( options: [:]                                   )
include { CENTRIFUGE                                          } from './modules/local/centrifuge'                  addParams( options: modules['centrifuge']                 )
include { KRAKEN2_DB_PREPARATION                              } from './modules/local/kraken2_db_preparation'      addParams( options: [:]                                   )
include { KRAKEN2                                             } from './modules/local/kraken2'                     addParams( options: modules['kraken2']                    )
include { KRONA_DB                                            } from './modules/local/krona_db'                    addParams( options: [:]                                   )
include { KRONA                                               } from './modules/local/krona'                       addParams( options: modules['krona']                      )
include { POOL_SINGLE_READS                                   } from './modules/local/pool_single_reads'           addParams( options: [:]                                   )
include { POOL_PAIRED_READS                                   } from './modules/local/pool_paired_reads'           addParams( options: [:]                                   )
include { POOL_SINGLE_READS as POOL_LONG_READS                } from './modules/local/pool_single_reads'           addParams( options: [:]                                   )
include { MEGAHIT                                             } from './modules/local/megahit'                     addParams( options: modules['megahit']                    )
include { SPADES                                              } from './modules/local/spades'                      addParams( options: modules['spades']                     )
include { SPADESHYBRID                                        } from './modules/local/spadeshybrid'                addParams( options: modules['spadeshybrid']               )
include { QUAST                                               } from './modules/local/quast'                       addParams( options: modules['quast']                      )
include { QUAST_BINS                                          } from './modules/local/quast_bins'                  addParams( options: modules['quast_bins']                 )
include { MERGE_QUAST_AND_BUSCO                               } from './modules/local/merge_quast_and_busco'       addParams( options: modules['merge_quast_and_busco']      )
include { CAT_DB                                              } from './modules/local/cat_db'                      addParams( options: [:]                                   )
include { CAT                                                 } from './modules/local/cat'                         addParams( options: modules['cat']                        )
include { MULTIQC                                             } from './modules/local/multiqc'                     addParams( options: multiqc_options                       )

// Local: Functions
include {
    hasExtension
} from './modules/local/functions'

// Local: Sub-workflows
include { METABAT2_BINNING    } from './subworkflows/local/metabat2_binning'      addParams( bowtie2_build_options: [:], bowtie2_align_options: modules['bowtie2_assembly_align'], metabat2_options: modules['metabat2']                                                   )
include { BUSCO_QC            } from './subworkflows/local/busco_qc'              addParams( busco_db_options: modules['busco_db_preparation'], busco_options: modules['busco'], busco_plot_options: modules['busco_plot'], busco_summary_options: modules['busco_summary'])

// nf-core/modules: Modules
include { FASTQC as FASTQC_RAW     } from './modules/nf-core/software/fastqc/main'              addParams( options: modules['fastqc_raw']            )
include { FASTQC as FASTQC_TRIMMED } from './modules/nf-core/software/fastqc/main'              addParams( options: modules['fastqc_trimmed']        )
include { FASTP                    } from './modules/nf-core/software/fastp/main'               addParams( options: modules['fastp']            )


////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config, params.phix_reference, params.host_fasta, params.centrifuge_db, params.kraken2_db, params.cat_db, params.lambda_reference, params.busco_reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.manifest) exit 1, "The parameter `--manifest` is deprecated.\n\tPlease use the `--input` parameter instead and mind the new format specifications."

// Check mandatory parameters
if (!params.input) exit 1, 'Input samplesheet not specified!'

hybrid = false
if(hasExtension(params.input, "csv")){
    // extracts read files from samplesheet CSV and distribute into channels
    Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() == 5) {
                    def id = row.sample
                    def group = row.group
                    def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : false
                    def sr2 = row.short_reads_2 ? file(row.short_reads_2, checkIfExists: true) : false
                    def lr = row.long_reads ? file(row.long_reads, checkIfExists: true) : false
                    if (lr) hybrid = true
                    // Check if given combination is valid
                    if (!sr1) exit 1, "Invalid input samplesheet: short_reads_1 can not be empty."
                    if (!sr2 && lr) exit 1, "Invalid input samplesheet: invalid combination of single-end short reads and long reads provided! SPAdes does not support single-end data and thus hybrid assembly cannot be performed."
                    if (!sr2 && !params.single_end) exit 1, "Invalid input samplesheet: single-end short reads provided, but command line parameter `--single_end` is false. Note that either only single-end or only paired-end reads must provided."
                    if (sr2 && params.single_end) exit 1, "Invalid input samplesheet: paired-end short reads provided, but command line parameter `--single_end` is true. Note that either only single-end or only paired-end reads must provided."
                    return [ id, group, sr1, sr2, lr ]
                } else {
                    exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                }
            }
        .set { ch_input_rows }
    // separate short and long reads
    ch_input_rows
        .map { id, group, sr1, sr2, lr ->
                    def meta = [:]
                    meta.id           = id
                    meta.group        = group
                    meta.single_end   = params.single_end
                    if (params.single_end) 
                        return [ meta, [ sr1] ]
                    else 
                        return [ meta, [ sr1, sr2 ] ]
            }
        .set { ch_raw_short_reads }
    ch_input_rows
        .map { id, group, sr1, sr2, lr ->
                    if (lr) {
                        def meta = [:]
                        meta.id           = id
                        meta.group        = group
                        return [ meta, lr ]
                    }
            }
        .set { ch_raw_long_reads }
} else {
    Channel
        .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .map { row -> 
                    def meta = [:]
                    meta.id           = row[0]
                    meta.group        = 0
                    meta.single_end   = params.single_end
                    return [ meta, row[1] ]
            }
        .set { ch_raw_short_reads }
    ch_input_rows = Channel.empty()
    ch_raw_long_reads = Channel.empty()
}

// Ensure sample IDs are unique
ch_input_rows
    .map { id, group, sr1, sr2, lr -> id }
    .toList()
    .map { ids -> if( ids.size() != ids.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample IDs!" } }

////////////////////////////////////////////////////
/* --     MORE PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check if binning mapping mode is valid
if (!['all','group','own'].contains(params.binning_map_mode))
    exit 1, "Invalid parameter '--binning_map_mode ${params.binning_map_mode}'. Valid values are 'all', 'group' or 'own'."
if (params.coassemble_group && params.binning_map_mode == 'own')
    exit 1, "Invalid combination of parameter '--binning_map_mode own' and parameter '--coassemble_group'. Select either 'all' or 'group' mapping mode when performing group-wise co-assembly."

// Check if specified cpus for SPAdes are available
if ( params.spades_fix_cpus && params.spades_fix_cpus > params.max_cpus )
    exit 1, "Invalid parameter '--spades_fix_cpus ${params.spades_fix_cpus}', max cpus are '${params.max_cpus}'."
if ( params.spadeshybrid_fix_cpus && params.spadeshybrid_fix_cpus > params.max_cpus )
    exit 1, "Invalid parameter '--spadeshybrid_fix_cpus ${params.spadeshybrid_fix_cpus}', max cpus are '${params.max_cpus}'."
// Check if settings concerning reproducibility of used tools are consistent and print warning if not
if (params.megahit_fix_cpu_1 || params.spades_fix_cpus || params.spadeshybrid_fix_cpus){
    if (!params.skip_spades && !params.spades_fix_cpus)
        log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes not. Consider using the parameter '--spades_fix_cpus'."
    if (hybrid && !params.skip_spadeshybrid && !params.spadeshybrid_fix_cpus)
        log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes hybrid not. Consider using the parameter '--spadeshybrid_fix_cpus'."
    if (!params.skip_megahit && !params.megahit_fix_cpu_1)
        log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but MEGAHIT not. Consider using the parameter '--megahit_fix_cpu_1'."
    if (!params.skip_binning && params.metabat_rng_seed == 0)
        log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but for MetaBAT2 a random seed is specified ('--metabat_rng_seed 0'). Consider specifying a positive seed instead."
}

// Check if SPAdes and single_end
if ( (!params.skip_spades || !params.skip_spadeshybrid) && params.single_end) {
    log.warn "metaSPAdes does not support single-end data. SPAdes will be skipped."
}

// Check if parameters for host contamination removal are valid and create channels
if ( params.host_fasta && params.host_genome) {
    exit 1, "Both host fasta reference and iGenomes genome are specififed to remove host contamination! Invalid combination, please specify either --host_fasta or --host_genome."
}
if ( hybrid && (params.host_fasta || params.host_genome) ) {
    log.warn "Host read removal is only applied to short reads. Long reads might be filtered indirectly by Filtlong, which is set to use read qualities estimated based on k-mer matches to the short, already filtered reads."
    if ( params.longreads_length_weight > 1 ) {
        log.warn "The parameter --longreads_length_weight is ${params.longreads_length_weight}, causing the read length being more important for long read filtering than the read quality. Set --longreads_length_weight to 1 in order to assign equal weights."
    }
}
if ( params.host_genome ) {
    // Check if host genome exists in the config file
    if ( !params.genomes.containsKey(params.host_genome) ) {
        exit 1, "The provided host genome '${params.host_genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    } else {
        host_fasta = params.genomes[params.host_genome].fasta ?: false
        if ( !host_fasta ) {
            exit 1, "No fasta file specified for the host genome ${params.host_genome}!"
        }
        Channel
            .value(file( "${host_fasta}" ))
            .set { ch_host_fasta }

        host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
        if ( !host_bowtie2index ) {
            exit 1, "No Bowtie 2 index file specified for the host genome ${params.host_genome}!"
        }
        Channel
            .value(file( "${host_bowtie2index}/*" ))
            .set { ch_host_bowtie2index }
    }
} else if ( params.host_fasta ) {
    Channel
        .value(file( "${params.host_fasta}" ))
        .set { ch_host_fasta }
} else {
    ch_host_fasta = Channel.empty()
}

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if(!params.skip_busco){
    Channel
        .value(file( "${params.busco_reference}" ))
        .set { ch_busco_db_file }
} else {
    ch_busco_db_file = Channel.empty()
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

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow {

    ch_software_versions = Channel.empty()

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    FASTQC_RAW (
        ch_raw_short_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))

    FASTP (
        ch_raw_short_reads
    )
    ch_short_reads = FASTP.out.reads
    ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))

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
        ch_software_versions = ch_software_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.version.first().ifEmpty(null))
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

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */
    NANOPLOT_RAW (
        ch_raw_long_reads
    )
    ch_software_versions = ch_software_versions.mix(NANOPLOT_RAW.out.version.first().ifEmpty(null))

    ch_long_reads = ch_raw_long_reads
    if (!params.skip_adapter_trimming) {
        PORECHOP (
            ch_raw_long_reads
        )
        ch_long_reads = PORECHOP.out.reads
        ch_software_versions = ch_software_versions.mix(PORECHOP.out.version.first().ifEmpty(null))
    }

    if (!params.keep_lambda) {
        NANOLYSE (
            ch_long_reads,
            ch_nanolyse_db
        )
        ch_long_reads = NANOLYSE.out.reads
        ch_software_versions = ch_software_versions.mix(NANOLYSE.out.version.first().ifEmpty(null))
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
    ch_software_versions = ch_software_versions.mix(FILTLONG.out.version.first().ifEmpty(null))

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
        CENTRIFUGE_DB_PREPARATION.out
    )
    ch_software_versions = ch_software_versions.mix(CENTRIFUGE.out.version.first().ifEmpty(null))

    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db_file
    )
    KRAKEN2 (
        ch_short_reads,
        KRAKEN2_DB_PREPARATION.out
    )
    ch_software_versions = ch_software_versions.mix(KRAKEN2.out.version.first().ifEmpty(null))

    if (( params.centrifuge_db || params.kraken2_db ) && !params.skip_krona){
        KRONA_DB ()
        KRONA (
            CENTRIFUGE.out.results_for_krona.mix(KRAKEN2.out.results_for_krona),
            KRONA_DB.out.collect()
        )
        ch_software_versions = ch_software_versions.mix(KRONA.out.version.first().ifEmpty(null))
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
        ch_software_versions = ch_software_versions.mix(MEGAHIT.out.version.first().ifEmpty(null))
    }

    // Co-assembly: pool reads for SPAdes
    if (params.coassemble_group) {
        // short reads
        if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)){
            if (params.single_end){
                POOL_SINGLE_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_SINGLE_READS.out
            } else {
                POOL_PAIRED_READS ( ch_short_reads_grouped )
                ch_short_reads_spades = POOL_PAIRED_READS.out
            }
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid){
            POOL_LONG_READS ( ch_long_reads_grouped )
            ch_long_reads_spades = POOL_LONG_READS.out
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
        ch_software_versions = ch_software_versions.mix(SPADES.out.version.first().ifEmpty(null))
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
        ch_software_versions = ch_software_versions.mix(SPADESHYBRID.out.version.first().ifEmpty(null))
    }

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast){
        QUAST ( ch_assemblies )
        ch_quast_multiqc = QUAST.out.qc
        ch_software_versions = ch_software_versions.mix(QUAST.out.version.first().ifEmpty(null))
    }

    /*
    ================================================================================
                                    Binning
    ================================================================================
    */

    ch_bowtie2_assembly_multiqc = Channel.empty()
    ch_busco_multiqc            = Channel.empty()
    if (!params.skip_binning){

        METABAT2_BINNING (
            ch_assemblies,
            ch_short_reads
        )
        ch_bowtie2_assembly_multiqc = METABAT2_BINNING.out.bowtie2_assembly_multiqc
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING.out.bowtie2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(METABAT2_BINNING.out.metabat2_version.first().ifEmpty(null))

        /*
        * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
        */
        BUSCO_QC (
            ch_busco_db_file,
            METABAT2_BINNING.out.bins.transpose()
        )
        ch_busco_multiqc = BUSCO_QC.out.multiqc
        ch_software_versions = ch_software_versions.mix(BUSCO_QC.out.version.first().ifEmpty(null))

        if (!params.skip_quast){
            QUAST_BINS ( METABAT2_BINNING.out.bins )
            ch_software_versions = ch_software_versions.mix(QUAST_BINS.out.version.first().ifEmpty(null))
            MERGE_QUAST_AND_BUSCO (
                QUAST_BINS.out.quast_bin_summaries.collect(),
                BUSCO_QC.out.summary
            )
        }

        /*
         * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
         */
        CAT_DB ( ch_cat_db_file )
        CAT ( 
            METABAT2_BINNING.out.bins,
            CAT_DB.out
        )
        ch_software_versions = ch_software_versions.mix(CAT.out.version.first().ifEmpty(null))
    }

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    /*
     * MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_custom_config.collect().ifEmpty([]),
            FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_removal_host_multiqc.collect{it[1]}.ifEmpty([]),
            ch_quast_multiqc.collect().ifEmpty([]),
            ch_bowtie2_assembly_multiqc.collect().ifEmpty([]),
            ch_busco_multiqc.collect().ifEmpty([])
        )
        multiqc_report       = MULTIQC.out.report.toList()
        //ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null)) // TODO ?
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
