/*
 * Binning with MetaBAT2 and MaxBin2
 */

params.mag_depths_options                           = [:]
params.mag_depths_plot_options                      = [:]
params.mag_depths_summary_options                   = [:]

include { METABAT2_METABAT2                     } from '../../modules/nf-core/modules/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../modules/nf-core/modules/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                               } from '../../modules/nf-core/modules/maxbin2/main'
include { GUNZIP                                } from '../../modules/nf-core/modules/gunzip/main'

include { CONVERT_DEPTHS                        } from '../../modules/local/convert_depths'
include { SPLIT_FASTA                           }  from '../../modules/local/split_fasta'
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

workflow BINNING {
    take:
    assemblies           // channel: [ val(meta), path(assembly), path(bams), path(bais) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:

    // TODO is scaffolds meant to go into here? These aren't being labelled correctly if so.
    ch_summarizedepth_input = assemblies
                                .map { meta, assembly, bams, bais ->
                                    [ meta, bams, bais ]
                                }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .set { ch_metabat_depths }

    ch_metabat2_input = assemblies
        .map { it ->
            it[0]['binner'] = 'MetaBAT2'

            [ it[0], it[1], it[2], it[3] ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { it ->
            [ it[0], it[1], it[4]]
        }

    // run binning, with convertion of metabat2 depth files to maxbin2 compatible if necessary
    // TODO consider replacing if statemetns with the new `when` modules args - requires latest modules though

    if ( !params.skip_metabat2 ) {
        METABAT2_METABAT2 ( ch_metabat2_input )
    }

    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( ch_metabat2_input )
            .map { it ->
                it[0]['binner'] = 'MaxBin2'

                [ it[0], it[1], it[2], it[3] ]
            }
            .set { ch_maxbin2_input }

        MAXBIN2 ( ch_maxbin2_input )
    }

    // split FASTQ
    if ( !params.skip_metabat2 & params.skip_maxbin2 ) {
        ch_input_splitfasta = METABAT2_METABAT2.out.unbinned
    } else if ( params.skip_metabat2 & !params.skip_maxbin2 ) {
        ch_input_splitfasta = MAXBIN2.out.unbinned_fasta
    } else {
        ch_input_splitfasta = METABAT2_METABAT2.out.unbinned.mix(MAXBIN2.out.unbinned_fasta)
    }

    SPLIT_FASTA ( ch_input_splitfasta )

    // decompress main bins (and large unbinned contigs) for MAG Depths,
    // first have to separate and re-group due to limitation of GUNZIP
    METABAT2_METABAT2.out.fasta.transpose().set { ch_metabat2_results_transposed }
    MAXBIN2.out.binned_fastas.transpose().set { ch_maxbin2_results_transposed }

    SPLIT_FASTA.out.unbinned.transpose().set { ch_split_fasta_results_transposed }

    ch_metabat2_results_transposed
        .mix( ch_maxbin2_results_transposed, ch_split_fasta_results_transposed )
        .set { ch_final_bins_for_gunzip }

    GUNZIP ( ch_final_bins_for_gunzip )
    GUNZIP.out.gunzip.groupTuple(by: 0).set{ ch_binning_results_gunzipped }

    // Compute bin depths for different samples (according to `binning_map_mode`)
    ch_depth_input = ch_binning_results_gunzipped
        .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0  )

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
    bins                                         = ch_metabat_results_gunzipped
    unbinned                                     = SPLIT_FASTA.out.unbinned
    tooshort                                     = METABAT2_METABAT2.out.tooshort
    lowdepth                                     = METABAT2_METABAT2.out.lowdepth
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    metabat2_version                             = METABAT2_METABAT2.out.versions
    metabat2_jgisummarizebamcontigdepths_version = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions
}
