include { MAG_DEPTHS                            } from '../../modules/local/mag_depths'
include { MAG_DEPTHS_PLOT                       } from '../../modules/local/mag_depths_plot'
include { MAG_DEPTHS_SUMMARY                    } from '../../modules/local/mag_depths_summary'

/*
 * Get number of columns in file (first line)
 */
def getColNo(filename) {
    lines  = file(filename).readLines()
    return lines[0].split('\t').size()
}

workflow DEPTHS {
    take:
    bins_unbins     //channel: val(meta), [ path(bins) ]
    depths          //channel: val(meta), path(depths)
    reads           //channel: val(meta), path(reads)

    main:
    ch_versions = Channel.empty()

    // Compute bin depths for different samples (according to `binning_map_mode`)
    // Create a new meta joining key first, but copy meta so that
    // we retain the information about binners and domain classification
    ch_depth_input = bins_unbins
        .map { meta, bins ->
            def meta_join = meta - meta.subMap('binner','domain')
            [ meta_join, meta, bins ]
        }
        .combine( depths, by: 0 )
        .map { meta_join, meta, bins, contig_depths_file ->
            [ meta, bins, contig_depths_file ]
        }
        .transpose()
        .groupTuple(by: [0,2])


    MAG_DEPTHS ( ch_depth_input )
    ch_versions = ch_versions.mix(MAG_DEPTHS.out.versions)

    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    // create file containing group information for all samples
    ch_sample_groups = reads
        .collectFile(name:'sample_groups.tsv'){ meta, reads -> meta.id + '\t' + meta.group + '\n' }

    // Filter MAG depth files: use only those for plotting that contain depths for > 2 samples
    ch_mag_depths_plot = MAG_DEPTHS.out.depths
        .map { meta, bin_depths_file ->
            if (getColNo(bin_depths_file) > 2) [ meta, bin_depths_file ]
        }

    MAG_DEPTHS_PLOT ( ch_mag_depths_plot, ch_sample_groups.collect() )
    MAG_DEPTHS_SUMMARY ( MAG_DEPTHS.out.depths.map{it[1]}.collect() )
    ch_versions = ch_versions.mix( MAG_DEPTHS_PLOT.out.versions )
    ch_versions = ch_versions.mix( MAG_DEPTHS_SUMMARY.out.versions )

    emit:
    depths_summary  = MAG_DEPTHS_SUMMARY.out.summary
    versions        = ch_versions
}
