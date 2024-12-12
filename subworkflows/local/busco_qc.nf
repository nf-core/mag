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
            BUSCO_UNTAR([[id: 'busco_db'], busco_db])

            if (busco_db.getSimpleName().contains('odb')) {
                busco_lineage = busco_db.getSimpleName()
            }
            busco_db = BUSCO_UNTAR.out.untar.map { it[1] }
        }
        else if (busco_db.isDirectory()) {
            if (busco_db.name.contains('odb')) {
                busco_lineage = busco_db.name
            }
        }
    }

    BUSCO_BUSCO( bins, 'genome', busco_lineage, busco_db, [] )

    // if (params.save_busco_db) {
    //     // publish files downloaded by Busco
    //     ch_downloads = BUSCO.out.busco_downloads.groupTuple().map { lin, downloads -> downloads[0] }.toSortedList().flatten()
    //     BUSCO_SAVE_DOWNLOAD(ch_downloads)
    // }

    COMBINE_BUSCO_TSV(BUSCO_BUSCO.out.batch_summary.map { it[1] }.collect())

    emit:
    summary    = COMBINE_BUSCO_TSV.out.combined
    multiqc    = BUSCO_BUSCO.out.short_summaries_txt.map { it[1] }.flatten()
    versions   = BUSCO_BUSCO.out.versions
}
