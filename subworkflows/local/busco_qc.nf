/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_DB_PREPARATION            } from '../../modules/local/busco_db_preparation'
include { BUSCO                           } from '../../modules/local/busco'
include { BUSCO_SAVE_DOWNLOAD             } from '../../modules/local/busco_save_download'
include { BUSCO_SUMMARY                   } from '../../modules/local/busco_summary'

workflow BUSCO_QC {
    take:
    busco_db   // channel: path
    bins       // channel: [ val(meta), path(bin) ]

    main:
    if ( busco_db.extension == 'gz' ) {
        // Expects to be tar.gz!
        BUSCO_DB_PREPARATION ( busco_db )

        ch_db_for_busco = BUSCO_DB_PREPARATION.out.db
                            .map{
                                meta, db ->
                                    def meta_new = [:]
                                    meta_new['id'] = meta
                                    meta_new['lineage'] = 'Y'
                                    [ meta_new, db ]
                            }
    } else if ( busco_db.isDirectory() ) {
        // Set meta to match expected channel cardinality for BUSCO
        ch_db_for_busco = Channel
                        .of(busco_db)
                        .map{
                            db ->
                                def meta = [:]
                                meta['id'] = db.toString().split('/').last()
                                if ("${meta['id'].toString().contains('odb10')}" == true) {
                                    meta['lineage'] = 'Y'
                                } else {
                                    meta['lineage'] = 'N'
                                }
                                [ meta, db ]
                        }
                        .collect()
    } else {
        error("Unsupported object given to --busco_db, database must be supplied as either a directory or a .tar.gz file!")
    }

    BUSCO (
        bins,
        ch_db_for_busco
    )

    if (params.save_busco_db){
        // publish files downloaded by Busco
        ch_downloads = BUSCO.out.busco_downloads.groupTuple().map{lin,downloads -> downloads[0]}.toSortedList().flatten()
        BUSCO_SAVE_DOWNLOAD ( ch_downloads )
    }

    BUSCO_SUMMARY (
        BUSCO.out.summary_domain.map{it[1]}.collect().ifEmpty([]),
        BUSCO.out.summary_specific.map{it[1]}.collect().ifEmpty([]),
        BUSCO.out.failed_bin.map{it[1]}.collect().ifEmpty([])
    )

    emit:
    summary     = BUSCO_SUMMARY.out.summary
    failed_bin  = BUSCO.out.failed_bin.map{it[1]}
    multiqc     = BUSCO.out.summary_domain.mix(BUSCO.out.summary_specific).map{it[1]}
    versions    = BUSCO.out.versions
}
