include { TIARA_TIARA                                                  } from '../../../modules/nf-core/tiara/tiara/main'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_TIARA } from '../../../modules/nf-core/dastool/fastatocontig2bin/main'
include { QSV_CAT as CONCAT_TIARA_TSV                                  } from '../../../modules/nf-core/qsv/cat/main'

include { TIARA_CLASSIFY                                               } from '../../../modules/local/tiara_classify/main'

workflow TIARA {
    take:
    ch_assemblies // [val(meta), path(fasta)]
    ch_in_bins    // [val(meta), path(fasta)]
    ch_in_unbins  // [val(meta), path(fasta)]

    main:
    ch_versions = channel.empty()

    ch_bins = ch_in_bins.map { meta, bin_list ->
        def meta_new = meta + [bin: 'bins']
        meta_new.bin = 'bins'
        [meta_new, bin_list]
    }

    ch_unbins = ch_in_unbins.map { meta, unbin_list ->
        def meta_new = meta + [bin: 'unbins']
        [meta_new, unbin_list]
    }

    ch_tiara_input = ch_bins.mix(ch_unbins)

    TIARA_TIARA(ch_assemblies)
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions)

    // Need contig2bin file for each bin group
    DASTOOL_FASTATOCONTIG2BIN_TIARA(ch_tiara_input, 'fa')
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_TIARA.out.versions)

    // Need to per-assembly Tiara classifications to their bins
    // Have to remove binner information from the meta map to do this
    ch_contigs_to_bin_tiara = DASTOOL_FASTATOCONTIG2BIN_TIARA.out.fastatocontig2bin
        .combine(ch_tiara_input, by: 0)
        .map { meta, contig2bin, bin_list ->
            def meta_join = meta - meta.subMap('binner', 'bin')
            [meta_join, meta, contig2bin, bin_list]
        }

    ch_tiara_classify_input = ch_contigs_to_bin_tiara
        .combine(TIARA_TIARA.out.classifications, by: 0)
        .map { _meta_join, meta, contig2bin, bin_list, classifications ->
            [meta, classifications, contig2bin, bin_list]
        }

    TIARA_CLASSIFY(ch_tiara_classify_input)
    ch_versions = ch_versions.mix(TIARA_CLASSIFY.out.versions)

    ch_eukarya_bins = TIARA_CLASSIFY.out.eukarya_bins.map { meta, bin_list ->
        def meta_new = meta + [domain: 'eukarya']
        [meta_new, bin_list]
    }

    ch_prokarya_bins = TIARA_CLASSIFY.out.prokarya_bins.map { meta, bin_list ->
        def meta_new = meta + [domain: 'prokarya']
        [meta_new, bin_list]
    }

    ch_bacteria_bins = TIARA_CLASSIFY.out.bacteria_bins.map { meta, bin_list ->
        def meta_new = meta + [domain: 'bacteria']
        [meta_new, bin_list]
    }

    ch_archaea_bins = TIARA_CLASSIFY.out.archaea_bins.map { meta, bin_list ->
        def meta_new = meta + [domain: 'archaea']
        [meta_new, bin_list]
    }

    ch_organelle_bins = TIARA_CLASSIFY.out.organelle_bins.map { meta, bin_list ->
        def meta_new = meta + [domain: 'organelle']
        [meta_new, bin_list]
    }

    ch_unknown_bins = TIARA_CLASSIFY.out.unknown_bins.map { meta, bin_list ->
        def meta_new = meta + [domain: 'unknown']
        [meta_new, bin_list]
    }

    ch_classified_bins_unbins = ch_eukarya_bins
        .mix(ch_prokarya_bins)
        .mix(ch_bacteria_bins)
        .mix(ch_archaea_bins)
        .mix(ch_organelle_bins)
        .mix(ch_unknown_bins)

    ch_classified_bins = ch_classified_bins_unbins
        .filter { meta, _bin_list ->
            meta.bin == "bins"
        }
        .map { meta, bin_list ->
            def meta_new = meta - meta.subMap('bin')
            [meta_new, bin_list]
        }

    ch_classified_unbins = ch_classified_bins_unbins
        .filter { meta, _bin_list ->
            meta.bin == "unbins"
        }
        .map { meta, bin_list ->
            def meta_new = meta - meta.subMap('bin')
            [meta_new, bin_list]
        }

    ch_bin_classifications = TIARA_CLASSIFY.out.bin_classifications
        .collect { _meta, classification ->
            classification
        }
        .map { bin_classifications ->
            [[id: 'tiara'], bin_classifications]
        }

    CONCAT_TIARA_TSV(ch_bin_classifications, 'rowskey', 'tsv', false)

    emit:
    classified_bins   = ch_classified_bins
    classified_unbins = ch_classified_unbins
    versions          = ch_versions
}
