/*
 * CheckM: Quantitative measures for the assessment of genome assembly
 */

include { CHECKM_QA        } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf/main'

workflow CHECKM_QC {
    take:
    bins       // channel: [ val(meta), path(bin) ]
    checkm_db  // TODO do proper input checks in mag.nf! Add to checkParams etc.,

    main:
    ch_versions = Channel.empty()

    ch_input_checkmdb = checkm_db ? checkm_db : []
    ch_bins_for_checkmlineagewf = bins
                                    .multiMap {
                                        meta, fa ->
                                            reads: [ meta, fa ]
                                            ext: fa.extension
                                    }

    CHECKM_LINEAGEWF ( ch_bins_for_checkmlineagewf.reads, ch_bins_for_checkmlineagewf.ext, checkm_db )
    ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)

    ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
        .join(CHECKM_LINEAGEWF.out.marker_file)
        .map{
            meta, dir, marker ->
            [ meta, dir, marker, []]
        }

    CHECKM_QA ( ch_checkmqa_input, [] )
    ch_versions = ch_versions.mix(CHECKM_QA.out.versions)

    // TODO Check output files published correctly

    emit:
    checkm_tsv = CHECKM_QA.out.output
    versions = ch_versions
}
