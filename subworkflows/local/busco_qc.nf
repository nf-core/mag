/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_SAVE_DOWNLOAD              } from '../../modules/local/busco_save_download'
include { BUSCO_BUSCO                      } from '../../modules/nf-core/busco/busco/main'
include { UNTAR as BUSCO_UNTAR             } from '../../modules/nf-core/untar'
include { COMBINE_TSV as COMBINE_BUSCO_TSV } from '../../modules/local/combine_tsv'

workflow BUSCO_QC {
    take:
    bins           // channel: [ val(meta), path(bin) ]
    busco_db       // channel: path
    busco_lineage  // channel: val

    main:
    if (!busco_db.isEmpty()) {
        if (busco_db.extension in ['gz', 'tgz']) {
            BUSCO_UNTAR(busco_db.map { db -> [[id: 'busco_db'], db] })

            busco_lineage = busco_db.getSimpleName()
            busco_db = BUSCO_UNTAR.out.untar.map { it[1] }

            // Expects to be tar.gz!
            // ch_db_for_busco = BUSCO_DB_PREPARATION(busco_db).db.map { meta, db ->
            //     def meta_new = [:]
            //     meta_new['id'] = meta
            //     meta_new['lineage'] = 'Y'
            //     [meta_new, db]
            // }
        }
        // else if (busco_db.isDirectory()) {
        //     // Set meta to match expected channel cardinality for BUSCO
        //     ch_db_for_busco = Channel
        //         .of(busco_db)
        //         .map { db ->
        //             def meta = [:]
        //             meta['id'] = db.getBaseName()
        //             if (meta['id'].contains('odb10') == true) {
        //                 meta['lineage'] = 'Y'
        //             }
        //             else {
        //                 meta['lineage'] = 'N'
        //             }
        //             [meta, db]
        //         }
        //         .collect()
        // }
        // BUSCO_BUSCO(bins, params.busco_lineage)
    }
    // else {
        // Set BUSCO database to empty to allow for --auto-lineage
    //     ch_db_for_busco = Channel
    //         .of([])
    //         .map { empty_db ->
    //             def meta = [:]
    //             meta['lineage'] = ''
    //             [meta, []]
    //         }
    //         .collect()
    // }

    BUSCO_BUSCO( bins, 'genome', busco_lineage, busco_db, [] )

    // if (params.save_busco_db) {
    //     // publish files downloaded by Busco
    //     ch_downloads = BUSCO.out.busco_downloads.groupTuple().map { lin, downloads -> downloads[0] }.toSortedList().flatten()
    //     BUSCO_SAVE_DOWNLOAD(ch_downloads)
    // }

    COMBINE_BUSCO_TSV(BUSCO_BUSCO.out.batch_summary.collect())

    emit:
    summary    = COMBINE_BUSCO_TSV.out.combined
    multiqc    = BUSCO_BUSCO.out.short_summaries_txt.map { it[1] }.flatten()
    versions   = BUSCO_BUSCO.out.versions
}
