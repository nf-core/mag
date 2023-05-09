/*
* Domain classification with Tiara
*/

include { TIARA } from '../../subworkflows/local/tiara'

workflow DOMAIN_CLASSIFICATION {
    take:
    assemblies // tuple val(meta), path(assembly)
    bins       // tuple val(meta), path( [ bins ] )
    unbins     // tuple val(meta), path( [ unbins ] )

    main:
    ch_versions = Channel.empty()

    if ( params.bin_domain_classification_tool == "tiara") {
        TIARA (assemblies, bins, unbins)
    }

    ch_classified_bins = TIARA.out.classified_bins
    ch_classified_unbins = TIARA.out.classified_unbins
    ch_versions = ch_versions.mix(TIARA.out.versions)

    emit:
    classified_bins   = ch_classified_bins
    classified_unbins = ch_classified_unbins
    versions          = ch_versions
}
