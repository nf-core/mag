/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN } from '../../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL                                        } from '../../../modules/nf-core/dastool/dastool/main.nf'

include { RENAME_PREDASTOOL                                      } from '../../../modules/local/dastool_rename_pre/main'
include { RENAME_POSTDASTOOL                                     } from '../../../modules/local/dastool_rename_post/main'

/*
 * Get number of columns in file (first line)
 */

workflow BINNING_REFINEMENT {
    take:
    ch_contigs_for_dastool // [val(meta), path(contigs)]
    ch_in_bins // [val(meta), path(bins)]

    main:
    ch_versions = channel.empty()

    // remove domain information, will add it back later
    // everything here is either unclassified or a prokaryote
    ch_bins = ch_in_bins
        .map { meta, bin_list ->
            def meta_new = meta - meta.subMap(['domain', 'refinement'])
            [meta_new, bin_list]
        }
        .groupTuple()
        .map { meta, bin_list ->
            [meta, bin_list.flatten()]
        }

    // prepare bins
    RENAME_PREDASTOOL(ch_bins)
    ch_bins_for_fastatocontig2bin = RENAME_PREDASTOOL.out.renamed_bins
        .transpose()
        .multiMap { meta, bin ->
            bins: [meta, bin]
            ext: bin.extension
        }
    ch_versions = ch_versions.mix(RENAME_PREDASTOOL.out.versions)

    // Generate DAS Tool auxilary files
    DASTOOL_FASTATOCONTIG2BIN(ch_bins_for_fastatocontig2bin.bins, ch_bins_for_fastatocontig2bin.ext)
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN.out.versions)

    // Run DAS Tool
    ch_fastatocontig2bin_for_dastool = DASTOOL_FASTATOCONTIG2BIN.out.fastatocontig2bin
        .map { meta, fastatocontig2bin ->
            def meta_new = meta - meta.subMap('binner')
            [meta_new, fastatocontig2bin]
        }
        .groupTuple(by: 0)

    // Note: do not `failOnMismatch` on join here, in some cases e.g. MAXBIN2 will fail if no bins, so cannot join!
    // Only want to join for DAS Tool on bins that 'exist'
    ch_input_for_dastool = ch_contigs_for_dastool
        .join(ch_fastatocontig2bin_for_dastool, by: 0)
        .map { meta, contigs, bins -> [meta, contigs, bins, []] }

    // Run DAS Tool
    DASTOOL_DASTOOL(ch_input_for_dastool, [])
    ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions)

    // Prepare bins for downstream analysis (separate from unbins, add 'binner' info and group)
    // use DAS Tool as 'binner' info allowing according grouping of refined bin sets,
    // while keeping information about original binning method in filenames and used binnames, e.g. "*-MaxBin2Refined-*.fa"
    // (alternatively one could think of adding, for example, meta.orig_binner, if this would simplify code)
    ch_dastool_bins_newmeta = DASTOOL_DASTOOL.out.bins
        .transpose()
        .map { meta, bin ->
            if (bin.name != "unbinned.fa") {
                def meta_new = meta + [binner: 'DASTool']
                [meta_new, bin]
            }
        }
        .groupTuple()
        .map { meta, bin_list ->
            def domain_class = params.bin_domain_classification ? 'prokarya' : 'unclassified'
            def meta_new = meta + [refinement: 'dastool_refined', domain: domain_class]
            [meta_new, bin_list]
        }

    ch_input_for_renamedastool = DASTOOL_DASTOOL.out.bins.map { meta, bin_list ->
        def domain_class = params.bin_domain_classification ? 'prokarya' : 'unclassified'
        def meta_new = meta + [refinement: 'dastool_refined', binner: 'DASTool', domain: domain_class]
        [meta_new, bin_list]
    }

    RENAME_POSTDASTOOL(ch_input_for_renamedastool)
    ch_versions = ch_versions.mix(RENAME_POSTDASTOOL.out.versions)

    refined_unbins = RENAME_POSTDASTOOL.out.refined_unbins.map { meta, bin_list ->
        def meta_new = meta + [refinement: 'dastool_refined_unbinned']
        [meta_new, bin_list]
    }

    emit:
    refined_bins   = ch_dastool_bins_newmeta
    refined_unbins = refined_unbins
    versions       = ch_versions
}
