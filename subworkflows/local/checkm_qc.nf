/*
 * CheckM/CheckM2: Quantitative measures for the assessment of genome assembly
 */

include { CHECKM_QA                         } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF                  } from '../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT                   } from '../../modules/nf-core/checkm2/predict/main'
include { COMBINE_TSV as COMBINE_CHECKM_TSV } from '../../modules/local/combine_tsv'

workflow CHECKM_QC {
    take:
    bins       // channel: [ val(meta), path(bin) ]
    checkm_db
    checkm2_db

    main:
    ch_versions = Channel.empty()

    if (params.binqc_tool == "checkm") {
        ch_bins_for_checkmlineagewf = bins.multiMap {
            meta, fa ->
                reads: [ meta, fa ]
                ext: fa.extension.unique().join("") // we set this in the pipeline to always `.fa` so this should be fine
        }

        CHECKM_LINEAGEWF(ch_bins_for_checkmlineagewf.reads, ch_bins_for_checkmlineagewf.ext, checkm_db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map{
                meta, dir, marker ->
                [ meta, dir, marker, []]
            }

        CHECKM_QA ( ch_checkmqa_input, [] )

        ch_versions = ch_versions.mix(CHECKM_QA.out.versions.first())

        COMBINE_CHECKM_TSV(CHECKM_QA.out.output.map{it[1]}.collect())
    }
    if (params.binqc_tool == "checkm2") {
        CHECKM2_PREDICT(bins, checkm2_db)

        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())

        COMBINE_CHECKM_TSV(CHECKM2_PREDICT.out.checkm2_tsv.map{it[1]}.collect())
    }

    emit:
    summary    = COMBINE_CHECKM_TSV.out.combined
    checkm_tsv = params.binqc_tool == "checkm" ? CHECKM_QA.out.output : []
    versions   = ch_versions
}
