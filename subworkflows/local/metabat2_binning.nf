/*
 * Binning with MetaBAT2
 */

params.mag_depths_options                           = [:]
params.mag_depths_plot_options                      = [:]
params.mag_depths_summary_options                   = [:]

include { METABAT2_METABAT2                  } from '../../modules/nf-core/modules/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../modules/nf-core/modules/metabat2/jgisummarizebamcontigdepths/main'
include { GUNZIP } from '../../modules/nf-core/modules/gunzip/main'

include { SPLIT_FASTQ }  from '../../modules/local/split_fastq'
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

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    ch_metabat2_input = assemblies
        .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0 )
        . map { it -> [ it[0], it[1], it[4] ] }

    METABAT2_METABAT2 ( ch_metabat2_input )

    METABAT2_METABAT2.out.fasta.dump(tag:"fasta_channel_metabat2_output")

    // split FASTQ
    SPLIT_FASTQ ( METABAT2_METABAT2.out.unbinned )

    // decompress main bins for downstream, have to separate and re-group due to limitation of GUNZIP
    METABAT2_METABAT2.out.fasta
        .transpose()
        .set { ch_metabat2_results_transposed }

    GUNZIP ( ch_metabat2_results_transposed )

    GUNZIP.out.gunzip.groupTuple(by: 0).set{ ch_metabat_results_gunzipped }

    // Compute bin depths for different samples (according to `binning_map_mode`)
    ch_depth_input = ch_metabat_results_gunzipped
        .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0  )
        .dump(tag:"pre_for_mag_depths")

    MAG_DEPTHS ( ch_depth_input )

    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    // create file containing group information for all samples
    ch_sample_groups = reads
        .collectFile(name:'sample_groups.tsv'){ meta, reads -> meta.id + '\t' + meta.group + '\n' }

    // filter MAG depth files: use only those for plotting that contain depths for > 2 samples
    ch_mag_depths_plot = MAG_DEPTHS.out.depths
        .map { meta, depth_file -> if (getColNo(depth_file) > 2) [meta, depth_file] }

    MAG_DEPTHS_PLOT ( ch_mag_depths_plot, ch_sample_groups.collect() )

    MAG_DEPTHS_SUMMARY ( MAG_DEPTHS.out.depths.map{it[1]}.collect() )

    emit:
    bins                                         = ch_metabat_results_gunzipped // TODO this would include discarded FASTAS! need to separate out somehow! Probably in metabat2 module
    unbinned                                     = SPLIT_FASTQ.out.unbinned
    tooshort                                     = METABAT2_METABAT2.out.tooshort
    lowdepth                                     = METABAT2_METABAT2.out.lowdepth
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    metabat2_version                             = METABAT2_METABAT2.out.versions
    metabat2_jgisummarizebamcontigdepths_version = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions
}
