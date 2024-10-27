/*
 * CheckM2: Assessing the quality of metagenome-derived genome bins using machine learning
 */

include { CHECKM2_PREDICT                       } from '../../modules/nf-core/checkm2/predict/main'
include { COMBINE_TSV as COMBINE_CHECKM2_TSV    } from '../../modules/local/combine_tsv'


workflow CHECKM2_QC {
    take:
    bins       // channel: [ val(meta), path(bin) ]
    checkm2_db

    main:
    ch_versions = Channel.empty()

    CHECKM2_PREDICT ( bins, checkm2_db )
    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())

    COMBINE_CHECKM2_TSV ( CHECKM2_PREDICT.out.checkm2_tsv.map{it[1]}.collect() )

    emit:
    summary     = COMBINE_CHECKM2_TSV.out.combined
    versions    = ch_versions
}
