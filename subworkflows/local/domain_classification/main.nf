/*
* Domain classification with Tiara
*/

include { TIARA } from '../../../subworkflows/local/tiara'

workflow DOMAIN_CLASSIFICATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_bins       // [val(meta), path(fasta)]
    ch_unbins     // [val(meta), path(fasta)]

    main:
    ch_versions = channel.empty()
    ch_tiara_summary = channel.empty()
    ch_publish_tiara = channel.empty()

    if (params.bin_domain_classification_tool == "tiara") {
        TIARA(ch_assemblies, ch_bins, ch_unbins)
        ch_versions = ch_versions.mix(TIARA.out.versions)
        ch_tiara_summary = TIARA.out.summary
        ch_publish_tiara = TIARA.out.publish_tiara
    }

    ch_classified_bins = TIARA.out.classified_bins
    ch_classified_unbins = TIARA.out.classified_unbins

    emit:
    classified_bins   = ch_classified_bins
    classified_unbins = ch_classified_unbins
    versions          = ch_versions
    tiara_summary     = ch_tiara_summary // [val(meta), path(tsv)]
    publish_tiara     = ch_publish_tiara // [val(meta), path(txt)]
}
