/*
* Domain classification with Tiara
*/

include { TIARA_TIARA                                                  } from '../../modules/nf-core/tiara/tiara/main'
include { TIARA_CLASSIFY                                               } from '../../modules/local/tiara_classify'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_TIARA } from '../../modules/nf-core/dastool/fastatocontig2bin/main'

workflow DOMAIN_CLASSIFICATION {
    take:
    assemblies // tuple val(meta), path(assembly)
    bins       // tuple val(meta), path( [ bins ])

    main:
    ch_versions = Channel.empty()

    TIARA_TIARA ( assemblies )
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions.first())
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

    ch_tiara_classify_input = ch_contigs_to_bin_tiara
        .combine( TIARA_TIARA.out.classifications, by: 0)
        .map { meta_join, meta, contig2bin, bins, classifications ->
            [ meta, classifications, contig2bin, bins ]
        }

    TIARA_CLASSIFY( ch_tiara_classify_input )
    ch_versions = ch_versions.mix(TIARA_CLASSIFY.out.versions.first())

    ch_eukarya_bins = TIARA_CLASSIFY.out.eukarya_bins
        .map { meta, bins ->
            def meta_new = meta.clone()
            meta_new.domain = 'eukarya'
            [meta_new, bins]
        }

    ch_prokarya_bins = TIARA_CLASSIFY.out.prokarya_bins
        .map { meta, bins ->
            def meta_new = meta.clone()
            meta_new.domain = 'prokarya'
            [meta_new, bins]
        }

    ch_bacteria_bins = TIARA_CLASSIFY.out.bacteria_bins
        .map { meta, bins ->
            def meta_new = meta.clone()
            meta_new.domain = 'bacteria'
            [meta_new, bins]
        }

    ch_archaea_bins = TIARA_CLASSIFY.out.archaea_bins
        .map { meta, bins ->
            def meta_new = meta.clone()
            meta_new.domain = 'archaea'
            [meta_new, bins]
        }

    ch_organelle_bins = TIARA_CLASSIFY.out.organelle_bins
        .map { meta, bins ->
            def meta_new = meta.clone()
            meta_new.domain = 'organelle'
            [meta_new, bins]
        }

    ch_unknown_bins = TIARA_CLASSIFY.out.unknown_bins
        .map { meta, bins ->
            def meta_new = meta.clone()
            meta_new.domain = 'unknown'
            [meta_new, bins]
        }

    ch_classified_bins = ch_eukarya_bins
        .mix(ch_prokarya_bins)
        .mix(ch_bacteria_bins)
        .mix(ch_archaea_bins)
        .mix(ch_organelle_bins)
        .mix(ch_unknown_bins)

    emit:
    classified_bins = ch_classified_bins
    versions        = ch_versions
}