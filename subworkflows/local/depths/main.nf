include { MAG_DEPTHS         } from '../../../modules/local/mag_depths/main'
include { MAG_DEPTHS_PLOT    } from '../../../modules/local/mag_depths_plot/main'
include { MAG_DEPTHS_SUMMARY } from '../../../modules/local/mag_depths_summary/main'

/*
 * Get number of columns in file (first line)
 */
def getColNo(filename) {
    def lines = file(filename).readLines()
    return lines[0].split('\t').size()
}

/*
 * Get number of rows in a file
 */
def getRowNo(filename) {
    def lines = file(filename).readLines()
    return lines.size()
}

workflow DEPTHS {
    take:
    ch_bins_unbins // [val(meta), path(fasta)]
    ch_depths      // [val(meta), path(depth)]
    ch_reads       // [val(meta), path(fastq)]

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

    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    // create file containing group information for all samples
    ch_sample_groups = ch_reads.collectFile(name: 'sample_groups.tsv') { meta, _sample_reads ->
        meta.id + '\t' + meta.group + '\n'
    }

    // Filter MAG depth files: use only those for plotting that contain depths for > 2 samples
    // as well as > 2 bins
    ch_mag_depths_plot = MAG_DEPTHS.out.depths.map { meta, bin_depths_file ->
        if (getColNo(bin_depths_file) > 2 && getRowNo(bin_depths_file) > 2) {
            [meta, bin_depths_file]
        }
    }

    MAG_DEPTHS_PLOT(ch_mag_depths_plot, ch_sample_groups.collect())
    ch_versions = ch_versions.mix(MAG_DEPTHS_PLOT.out.versions)

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
