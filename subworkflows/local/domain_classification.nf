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
    TIARA ( assemblies )
    // Need contig2bin file for each bin group
    DASTOOL_FASTATOCONTIG2BIN_TIARA ( bins , 'fa')

    // Need to per-assembly Tiara classifications to their bins
    // Have to remove binner information from the meta map to do this
    ch_contigs_to_bin_tiara = DASTOOL_FASTATOCONTIG2BIN_TIARA.out.fastatocontig2bin
        .map { meta, bins ->
            def meta_join = meta.clone()
            meta_join.remove('binner')
            [ meta_join, meta, bins ]
        }

    ch_tiara_classify_input = TIARA.out.classifications
        .combine(ch_contigs_to_bin_tiara, by: 0)
        .map { meta_join, classifications, meta, bins ->
            [ meta, classifications, bins ]
        }

    TIARA_CLASSIFY( ch_tiara_classify_input )

}
