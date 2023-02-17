/*
* Domain classification with Tiara
*/

include { TIARA                                                        } from '../../modules/local/tiara'
include { TIARA_CLASSIFY                                               } from '../../modules/local/tiara_classify'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_TIARA } from '../../modules/nf-core/dastool/fastatocontig2bin/main'

workflow DOMAIN_CLASSIFICATION {
    take:
    assemblies // tuple val(meta), path(assembly)
    bins       // tuple val(meta), path(bins)

    main:
    ch_versions = Channel.empty()

    TIARA ( assemblies )
    ch_versions = ch_versions.mix(TIARA.out.versions.first())
    // Need contig2bin file for each bin group
    DASTOOL_FASTATOCONTIG2BIN_TIARA ( bins , 'fa')
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_TIARA.out.versions.first())

    // Need to per-assembly Tiara classifications to their bins
    // Have to remove binner information from the meta map to do this
    ch_contigs_to_bin_tiara = DASTOOL_FASTATOCONTIG2BIN_TIARA.out.fastatocontig2bin
        .combine(bins, by: 0)
        .map { meta, contig2bin, bins ->
            def meta_join = meta.clone()
            meta_join.remove('binner')
            [ meta_join, meta, contig2bin, bins ]
        }

    ch_tiara_classify_input = TIARA.out.classifications
        .combine(ch_contigs_to_bin_tiara, by: 0)
        .map { meta_join, classifications, meta, contig2bin, bins ->
            [ meta, classifications, contig2bin, bins ]
        }

    TIARA_CLASSIFY( ch_tiara_classify_input )
    ch_versions = ch_versions.mix(TIARA_CLASSIFY.out.versions.first())

    emit:
    eukarya_bins   = TIARA_CLASSIFY.out.eukarya_bins
    prokarya_bins  = TIARA_CLASSIFY.out.prokarya_bins
    bacteria_bins  = TIARA_CLASSIFY.out.bacteria_bins
    archaea_bins   = TIARA_CLASSIFY.out.archaea_bins
    organelle_bins = TIARA_CLASSIFY.out.organelle_bins
    unknown_bins   = TIARA_CLASSIFY.out.unknown_bins
    versions       = ch_versions
}
