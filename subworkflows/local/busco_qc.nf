/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_DB_PREPARATION            } from '../../modules/local/busco_db_preparation'
include { BUSCO                           } from '../../modules/local/busco'
include { BUSCO_SAVE_DOWNLOAD             } from '../../modules/local/busco_save_download'
include { BUSCO_SUMMARY                   } from '../../modules/local/busco_summary'

workflow BUSCO_QC {
    take:
    busco_db_file           // channel: path
    busco_download_folder   // channel: path
    bins                    // channel: [ val(meta), path(bin) ]

    main:
    if (params.busco_reference){
        BUSCO_DB_PREPARATION ( busco_db_file )
        ch_busco_db = BUSCO_DB_PREPARATION.out.db
    } else {
        ch_busco_db = Channel.empty()
    }
    BUSCO (
        bins,
        ch_busco_db.collect().ifEmpty([]),
        busco_download_folder.collect().ifEmpty([])
    )
    if (params.save_busco_reference){
        // publish files downloaded by Busco
        ch_downloads = BUSCO.out.busco_downloads.groupTuple().map{lin,downloads -> downloads[0]}.toSortedList().flatten()
        BUSCO_SAVE_DOWNLOAD ( ch_downloads )
    }

    busco_summary_domain = BUSCO.out.summary_domain.collect()
    busco_summary_specific = BUSCO.out.summary_specific.collect()
    busco_failed_bin = BUSCO.out.failed_bin.collect()

    busco_summary_domain.dump(tag: 'busco_summary_domain', pretty: true)
    busco_summary_specific.dump(tag: 'busco_summary_specific', pretty: true)
    busco_failed_bin.dump(tag: 'busco_failed_bin', pretty: true)

    BUSCO_SUMMARY (
        BUSCO.out.summary_domain.map{it[1]}.unique().collect().ifEmpty([]),
        BUSCO.out.summary_specific.map{it[1]}.unique().collect().ifEmpty([]),
        BUSCO.out.failed_bin.map{it[1]}.unique().collect().ifEmpty([])
    )

    emit:
    summary     = BUSCO_SUMMARY.out.summary
    failed_bin  = BUSCO.out.failed_bin.map{it[1]}
    multiqc     = BUSCO.out.summary_domain.mix(BUSCO.out.summary_specific).map{it[1]}
    versions    = BUSCO.out.versions
}
