/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

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

    COMBINE_BUSCO_TSV(BUSCO_BUSCO.out.batch_summary.map { it[1] }.collect())

    versions = BUSCO_BUSCO.out.versions.first().mix(COMBINE_BUSCO_TSV.out.versions)

    emit:
    summary    = COMBINE_BUSCO_TSV.out.combined
    multiqc    = BUSCO_BUSCO.out.short_summaries_txt.map { it[1] }.flatten()
    versions   = versions
}
