include { ARIA2 as ARIA2_UNTAR } from '../../modules/nf-core/aria2/main'
include { BUSCO_QC             } from '../../subworkflows/local/busco_qc'
include { CHECKM_QC            } from '../../subworkflows/local/checkm_qc'
include { GUNC_QC              } from '../../subworkflows/local/gunc_qc'
include { QUAST_BINS           } from '../../modules/local/quast_bins'
include { QUAST_BINS_SUMMARY   } from '../../modules/local/quast_bins_summary'

workflow BIN_QC {
    take:
    bins_unbins

    main:
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

    // Get checkM database if not supplied
    if ( !params.skip_binqc && params.binqc_tool == 'checkm' && !params.checkm_db ) {
        ARIA2_UNTAR (params.checkm_download_url)
        ch_checkm_db = ARIA2_UNTAR.out.downloaded_file
    }

    if (params.gunc_db) {
    ch_gunc_db = file(params.gunc_db, checkIfExists: true)
    } else {
        ch_gunc_db = Channel.empty()
    }

    ch_versions = Channel.empty()

    ch_busco_multiqc = Channel.empty()
    ch_busco_summary = Channel.empty()
    ch_checkm_summary = Channel.empty()
    ch_quast_bins_summary = Channel.empty()

    bins_unbins_transposed = bins_unbins.transpose()

    if (!params.skip_binqc && params.binqc_tool == 'busco'){
        /*
        * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
        */
        BUSCO_QC (
            ch_busco_db_file,
            ch_busco_download_folder,
            bins_unbins_transposed
        )
        ch_busco_summary = BUSCO_QC.out.summary
        ch_busco_multiqc = BUSCO_QC.out.multiqc
        ch_versions = ch_versions.mix(BUSCO_QC.out.versions.first())
    }

    if (!params.skip_binqc && params.binqc_tool == 'checkm'){
        /*
        * CheckM subworkflow: Quantitative measures for the assessment of genome assembly
        */
        CHECKM_QC (
            bins_unbins_transposed.groupTuple(),
            ch_checkm_db
        )
        ch_checkm_summary = CHECKM_QC.out.summary

        // TODO custom output parsing? Add to MultiQC?
        ch_versions = ch_versions.mix(CHECKM_QC.out.versions)

    }

    if ( params.run_gunc && params.binqc_tool == 'checkm' ) {
        GUNC_QC ( bins_unbins_transposed, ch_gunc_db, CHECKM_QC.out.checkm_tsv )
        ch_versions = ch_versions.mix( GUNC_QC.out.versions )
    } else if ( params.run_gunc ) {
        GUNC_QC ( bins_unbins_transposed, ch_gunc_db, [] )
        ch_versions = ch_versions.mix( GUNC_QC.out.versions )
    }

    if (!params.skip_quast){
        ch_input_for_quast_bins = bins_unbins
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

    emit:
    busco_summary      = ch_busco_summary
    busco_multiqc      = ch_busco_multiqc
    busco_failed_bins  = BUSCO_QC.out.failed_bin
    checkm_summary     = ch_checkm_summary
    quast_bins_summary = ch_quast_bins_summary
    versions           = ch_versions
}
