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
    if ( !busco_db.isEmpty() ) {
        if ( busco_db.extension in ['gz', 'tgz'] ) {
            // Expects to be tar.gz!
            ch_db_for_busco = BUSCO_DB_PREPARATION ( busco_db ).db
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
                                        meta['id'] = db.getBaseName()
                                        if ( meta['id'].contains('odb10') == true ) {
                                            meta['lineage'] = 'Y'
                                        } else {
                                            meta['lineage'] = 'N'
                                        }
                                        [ meta, db ]
                                }
                                .collect()
        }
    } else {
        // Set BUSCO database to empty to allow for --auto-lineage
        ch_db_for_busco = Channel
                            .of([])
                            .map{
                                empty_db ->
                                    def meta = [:]
                                    meta['lineage'] = ''
                                    [ meta, [] ]
                            }
                            .collect()
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

    busco_summary_domain = BUSCO.out.summary_domain.collect()
    busco_summary_specific = BUSCO.out.summary_specific.collect()
    busco_failed_bin = BUSCO.out.failed_bin.collect()

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
