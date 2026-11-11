include { MAG_DEPTHS         } from '../../../modules/local/mag_depths/main'
include { MAG_DEPTHS_SUMMARY } from '../../../modules/local/mag_depths_summary/main'

workflow DEPTHS {
    take:
    ch_bins_unbins // [val(meta), path(fasta)]
    ch_depths      // [val(meta), path(depth)]

    main:
    ch_versions = channel.empty()

    // Compute bin depths for different samples (according to `binning_map_mode`)
    // Create a new meta combine key first, but copy meta so that
    // we retain the information about binners and domain classification
    ch_depth_input = ch_bins_unbins
        .map { meta, bins ->
            def meta_combine = meta - meta.subMap('binner', 'domain', 'refinement')
            [meta_combine, meta, bins]
        }
        .groupTuple()
        .combine(ch_depths, by: 0)
        .transpose()
        .map { _meta_combine, meta, bins, depth ->
            def meta_new = meta - meta.subMap('domain', 'refinement')
            [meta_new, bins, depth]
        }
        .groupTuple(by: [0, 2])
        .map { meta, bins, depth ->
            [meta, bins.unique().flatten(), depth]
        }

    MAG_DEPTHS(ch_depth_input)
    ch_versions = ch_versions.mix(MAG_DEPTHS.out.versions)

    //Depth files that are coming from bins and failed binning refinement are concatenated per meta
    ch_mag_depth_out = MAG_DEPTHS.out.depths.collectFile(keepHeader: true) { meta, depth ->
        [meta.id, depth]
    }

    MAG_DEPTHS_SUMMARY(ch_mag_depth_out.collect())
    ch_versions = ch_versions.mix(MAG_DEPTHS_SUMMARY.out.versions)

    emit:
    depths_summary = MAG_DEPTHS_SUMMARY.out.summary
    versions       = ch_versions
}
