/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_METABAT2 } from '../../modules/nf-core/modules/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_MAXBIN2  } from '../../modules/nf-core/modules/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL                                                 } from '../../modules/nf-core/modules/dastool/dastool/main.nf'
include { RENAME_PREDASTOOL                                               } from '../../modules/local/rename_predastool'
include { RENAME_POSTDASTOOL                                              } from '../../modules/local/rename_postdastool'
include { MAG_DEPTHS as MAG_DEPTHS_REFINED                                } from '../../modules/local/mag_depths'
include { MAG_DEPTHS_PLOT as MAG_DEPTHS_PLOT_REFINED                      } from '../../modules/local/mag_depths_plot'
include { MAG_DEPTHS_SUMMARY as MAG_DEPTHS_SUMMARY_REFINED                } from '../../modules/local/mag_depths_summary'

/*
 * Get number of columns in file (first line)
 */
def getColNo(filename) {
    lines  = file(filename).readLines()
    return lines[0].split('\t').size()
}

workflow BINNING_REFINEMENT {
    take:
    contigs        //
    bins           // channel: [ val(meta), path(bins) ]
    depths
    reads

    main:
    ch_versions = Channel.empty()

    // Drop unnecessary files
    ch_contigs_for_dastool = contigs
                                .map {
                                    meta, assembly, bams, bais ->
                                        def meta_new = meta.clone()
                                        [ meta_new, assembly ]
                                }

    ch_bins_for_fastatocontig2bin = RENAME_PREDASTOOL(bins).renamed_bins
                                        .branch {
                                            metabat2: it[0]['binner'] == 'MetaBAT2'
                                            maxbin2:  it[0]['binner'] == 'MaxBin2'
                                        }

    // Generate DASTool auxilary files
    DASTOOL_FASTATOCONTIG2BIN_METABAT2 ( ch_bins_for_fastatocontig2bin.metabat2, "fa")
    // MaxBin2 bin extension was changed to 'fa' as well in RENAME_PREDASTOOL
    DASTOOL_FASTATOCONTIG2BIN_MAXBIN2 ( ch_bins_for_fastatocontig2bin.maxbin2, "fa")

    // Run DASTOOL
    ch_fastatocontig2bin_for_dastool = Channel.empty()
    ch_fastatocontig2bin_for_dastool = ch_fastatocontig2bin_for_dastool
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_METABAT2.out.fastatocontig2bin)
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_MAXBIN2.out.fastatocontig2bin)
                                    .map {
                                        meta, fastatocontig2bin ->
                                            def meta_new = meta.clone()
                                            meta_new.remove('binner')
                                            [ meta_new, fastatocontig2bin ]
                                    }
                                    .groupTuple(by: 0)

    ch_input_for_dastool = ch_contigs_for_dastool.join(ch_fastatocontig2bin_for_dastool, by: 0)

    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_METABAT2.out.versions.first())
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_MAXBIN2.out.versions.first())

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
                    def meta_new = meta.clone()
                    meta_new['binner'] = 'DASTool'
                    [ meta_new, bin ]
                }
            }
        .groupTuple()

    ch_input_for_renamedastool = DASTOOL_DASTOOL.out.bins
        .map {
            meta, bins ->
                def meta_new = meta.clone()
                meta_new['binner'] = 'DASTool'
                [ meta_new, bins ]
            }

    RENAME_POSTDASTOOL ( ch_input_for_renamedastool )

    // We have to strip the meta to be able to combine with the original
    // depths file to run MAG_DEPTH
    ch_input_for_magdepth = ch_dastool_bins_newmeta
        .mix( RENAME_POSTDASTOOL.out.refined_unbins )
        .map {
                meta, refinedbins ->
                def meta_new = meta.clone()
                meta_new.remove('binner')
                [ meta_new, refinedbins ]
        }
        .transpose()
        .groupTuple (by: 0 )
        .join( depths, by: 0 )

    MAG_DEPTHS_REFINED ( ch_input_for_magdepth )

    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    // create file containing group information for all samples
    ch_sample_groups = reads
        .collectFile(name:'sample_groups.tsv'){ meta, reads -> meta.id + '\t' + meta.group + '\n' }

    // Transpose and add 'binner' meta information again for plotting
    // filter MAG depth files: use only those for plotting that contain depths for > 2 samples
    ch_mag_depths_plot_refined = MAG_DEPTHS_REFINED.out.depths
        .transpose()
        .map { meta, depth_file ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'DASTool'
            if (getColNo(depth_file) > 2) [meta_new, depth_file]
        }

    MAG_DEPTHS_PLOT_REFINED ( ch_mag_depths_plot_refined, ch_sample_groups.collect() )
    MAG_DEPTHS_SUMMARY_REFINED ( MAG_DEPTHS_REFINED.out.depths.map{it[1]}.collect() )

    emit:
    refined_bins                = ch_dastool_bins_newmeta
    refined_unbins              = RENAME_POSTDASTOOL.out.refined_unbins
    refined_depths              = MAG_DEPTHS_REFINED.out.depths
    refined_depths_summary      = MAG_DEPTHS_SUMMARY_REFINED.out.summary
    versions                    = ch_versions
}
