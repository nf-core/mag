/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_METABAT2 } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_MAXBIN2  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_CONCOCT  } from '../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL                                                 } from '../../modules/nf-core/dastool/dastool/main.nf'
include { RENAME_PREDASTOOL                                               } from '../../modules/local/rename_predastool'
include { RENAME_POSTDASTOOL                                              } from '../../modules/local/rename_postdastool'

/*
 * Get number of columns in file (first line)
 */

workflow BINNING_REFINEMENT {
    take:
    ch_contigs_for_dastool // channel: [ val(meta), path(contigs) ]
    bins           // channel: [ val(meta), path(bins) ]

    main:
    ch_versions = Channel.empty()

    // remove domain information, will add it back later
    // everything here is either unclassified or a prokaryote
    ch_bins = bins
        .map { meta, bins ->
            def meta_new = meta - meta.subMap(['domain','refinement'])
            [meta_new, bins]
        }
        .groupTuple()
        .map {
            meta, bins -> [meta, bins.flatten()]
        }

    // prepare bins
    ch_bins_for_fastatocontig2bin = RENAME_PREDASTOOL(ch_bins).renamed_bins
                                        .branch {
                                            metabat2: it[0]['binner'] == 'MetaBAT2'
                                            maxbin2:  it[0]['binner'] == 'MaxBin2'
                                            concoct:  it[0]['binner'] == 'CONCOCT'
                                        }

    // Generate DASTool auxilary files
    DASTOOL_FASTATOCONTIG2BIN_METABAT2 ( ch_bins_for_fastatocontig2bin.metabat2, "fa")
    // MaxBin2 bin extension was changed to 'fa' as well in RENAME_PREDASTOOL
    DASTOOL_FASTATOCONTIG2BIN_MAXBIN2 ( ch_bins_for_fastatocontig2bin.maxbin2, "fa")
    DASTOOL_FASTATOCONTIG2BIN_CONCOCT ( ch_bins_for_fastatocontig2bin.concoct, "fa")

    // Run DASTOOL
    ch_fastatocontig2bin_for_dastool = Channel.empty()
    ch_fastatocontig2bin_for_dastool = ch_fastatocontig2bin_for_dastool
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_METABAT2.out.fastatocontig2bin)
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_MAXBIN2.out.fastatocontig2bin)
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_CONCOCT.out.fastatocontig2bin)
                                    .map {
                                        meta, fastatocontig2bin ->
                                            def meta_new = meta - meta.subMap('binner')
                                            [ meta_new, fastatocontig2bin ]
                                    }
                                    .groupTuple(by: 0)

    // Note: do not `failOnMismatch` on join here, in some cases e.g. MAXBIN2 will fail if no bins, so cannot join!
    // Only want to join for DAS_Tool on bins that 'exist'

    ch_input_for_dastool = ch_contigs_for_dastool.join(ch_fastatocontig2bin_for_dastool, by: 0)

    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_METABAT2.out.versions.first())
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_MAXBIN2.out.versions.first())
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_CONCOCT.out.versions.first())

    // Run DAStool
    DASTOOL_DASTOOL(ch_input_for_dastool, [], [])
    ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions.first())

    // Prepare bins for downstream analysis (separate from unbins, add 'binner' info and group)
    // use DASTool as 'binner' info allowing according grouping of refined bin sets,
    // while keeping information about original binning method in filenames and used binnames, e.g. "*-MaxBin2Refined-*.fa"
    // (alternatively one could think of adding, for example, meta.orig_binner, if this would simplify code)
    ch_dastool_bins_newmeta = DASTOOL_DASTOOL.out.bins.transpose()
        .map {
            meta, bin ->
                if (bin.name != "unbinned.fa") {
                    def meta_new = meta + [binner: 'DASTool']
                    [ meta_new, bin ]
                }
            }
        .groupTuple()
        .map {
            meta, bins ->
                def domain_class = params.bin_domain_classification ? 'prokarya' : 'unclassified'
                def meta_new = meta + [refinement: 'dastool_refined', domain: domain_class]
                [ meta_new, bins ]
            }

    ch_input_for_renamedastool = DASTOOL_DASTOOL.out.bins
        .map {
            meta, bins ->
                def domain_class = params.bin_domain_classification ? 'prokarya' : 'unclassified'
                def meta_new = meta + [refinement: 'dastool_refined', binner: 'DASTool', domain: domain_class]
                [ meta_new, bins ]
            }

    RENAME_POSTDASTOOL ( ch_input_for_renamedastool )

    refined_unbins = RENAME_POSTDASTOOL.out.refined_unbins
        .map {
            meta, bins ->
                def meta_new = meta + [refinement: 'dastool_refined_unbinned']
                [meta_new, bins]
        }

    emit:
    refined_bins                = ch_dastool_bins_newmeta
    refined_unbins              = refined_unbins
    versions                    = ch_versions
}
