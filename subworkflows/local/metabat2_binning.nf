/*
 * Binning with MetaBAT2
 */

params.metabat2_jgisummarizebamcontigdepths_options = [:]
params.metabat2_options                             = [:]
params.mag_depths_options                           = [:]
params.mag_depths_plot_options                      = [:]
params.mag_depths_summary_options                   = [:]

include { METABAT2                  } from '../../modules/local/metabat2'                 addParams( options: params.metabat2_options           ) // local

include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../modules/nf-core/modules/metabat2/jgisummarizebamcontigdepths/main' addParams( options: [ suffix: { meta.assembler } ] )
include { MAG_DEPTHS                            } from '../../modules/local/mag_depths'               addParams( options: params.mag_depths_options         )
include { MAG_DEPTHS_PLOT                       } from '../../modules/local/mag_depths_plot'          addParams( options: params.mag_depths_plot_options    )
include { MAG_DEPTHS_SUMMARY                    } from '../../modules/local/mag_depths_summary'       addParams( options: params.mag_depths_summary_options )

/*
 * Get number of columns in file (first line)
 */
def getColNo(filename) {
    lines  = file(filename).readLines()
    return lines[0].split('\t').size()
}

workflow METABAT2_BINNING {
    take:
    assemblies           // channel: [ val(meta), path(assembly), path(bams), path(bais) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:

    ch_summarizedepth_input = assemblies
                                .map { meta, assembly, bams, bais -> [ meta, bams, bais ] }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS (
        ch_summarizedepth_input
    )

    // TODO re-merge with FASTA assemblies
    ch_metabat2_input = assemblies
        .join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth)

    METABAT2 ( ch_metabat2_input )

    // Compute bin depths for different samples (according to `binning_map_mode`)
    MAG_DEPTHS (
        METABAT2.out.bins,
        METABAT2.out.depths
    )
    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    // create file containing group information for all samples
    ch_sample_groups = reads
        .collectFile(name:'sample_groups.tsv'){ meta, reads -> meta.id + '\t' + meta.group + '\n' }

    // filter MAG depth files: use only those for plotting that contain depths for > 2 samples
    ch_mag_depths_plot = MAG_DEPTHS.out.depths
        .map { meta, depth_file -> if (getColNo(depth_file) > 2) [meta, depth_file] }

    MAG_DEPTHS_PLOT (
        ch_mag_depths_plot,
        ch_sample_groups.collect()
    )

    MAG_DEPTHS_SUMMARY ( MAG_DEPTHS.out.depths.map{it[1]}.collect() )

    emit:
    bins                                         = METABAT2.out.bins
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    metabat2_version                             = METABAT2.out.version
    metabat2_jgisummarizebamcontigdepths_version = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions
}
