/*
* Domain classification with Tiara
*/

include { TIARA } from '../../../subworkflows/local/tiara'

workflow DOMAIN_CLASSIFICATION {
    take:
    ch_assemblies // tuple val(meta), path(assembly)
    ch_bins       // tuple val(meta), path( [ bins ] )
    ch_unbins     // tuple val(meta), path( [ unbins ] )

    main:
    ch_versions = Channel.empty()

    if (params.bin_domain_classification_tool == "tiara") {
        TIARA(ch_assemblies, ch_bins, ch_unbins)
    }

    ch_classified_bins = TIARA.out.classified_bins
    ch_classified_unbins = TIARA.out.classified_unbins
    ch_versions = ch_versions.mix(TIARA.out.versions)

    emit:
    classified_bins   = ch_classified_bins
    classified_unbins = ch_classified_unbins
    versions          = ch_versions
}
