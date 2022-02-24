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
include { SPLIT_FASTA                           } from '../../modules/local/split_fasta'
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

    // TODO are scaffolds meant to go into here? These aren't being labelled
    // correctly if so...

    // generate coverage depths for each contig
    ch_summarizedepth_input = assemblies
                                .map { meta, assembly, bams, bais ->
                                        def meta_new = meta.clone()
                                    [ meta_new, bams, bais ]
                                }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, depths ]
        }
        .set { ch_metabat_depths }

    // combine depths back with assemblies
    ch_metabat2_input = assemblies
        .map { meta, contigs, reads, indicies ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, contigs, reads, indicies ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, contigs, reads, indicies, depths ->
            [ meta, contigs, depths ]
        }

    // conver metabat2 depth files to maxbin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( ch_metabat2_input )
        CONVERT_DEPTHS.out.output
            .map { meta, contigs, reads, depth ->
                    def meta_new = meta.clone()
                    meta_new['binner'] = 'MaxBin2'

                [ meta_new, contigs, reads, depth ]
            }
            .set { ch_maxbin2_input }
    }

    // run binning
    if ( !params.skip_metabat2 ) {
        METABAT2_METABAT2 ( ch_metabat2_input )
    }

    if ( !params.skip_maxbin2 ) {
        MAXBIN2 ( ch_maxbin2_input )
    }

    // split fastq files, depending
    if ( !params.skip_metabat2 & params.skip_maxbin2 ) {
        ch_input_splitfasta = METABAT2_METABAT2.out.unbinned
    } else if ( params.skip_metabat2 & !params.skip_maxbin2 ) {
        ch_input_splitfasta = MAXBIN2.out.unbinned_fasta
    } else {
        ch_input_splitfasta = METABAT2_METABAT2.out.unbinned.mix(MAXBIN2.out.unbinned_fasta)
    }

    SPLIT_FASTA ( ch_input_splitfasta )

    // decompress main bins (and large unbinned contigs from SPLIT_FASTA) for
    // MAG Depths, first have to separate and re-group due to limitation of
    // GUNZIP module
    METABAT2_METABAT2.out.fasta.transpose().set { ch_metabat2_results_transposed }
    MAXBIN2.out.binned_fastas.transpose().set { ch_maxbin2_results_transposed }
    SPLIT_FASTA.out.unbinned.transpose().set { ch_split_fasta_results_transposed }

    ch_metabat2_results_transposed
        .mix( ch_maxbin2_results_transposed, ch_split_fasta_results_transposed )
        .set { ch_final_bins_for_gunzip }

    GUNZIP ( ch_final_bins_for_gunzip )
    GUNZIP.out.gunzip
        .set{ ch_binning_results_gunzipped }

    // Compute bin depths for different samples (according to `binning_map_mode`)
    // Have to remove binner meta for grouping to mix back with original depth
    // files, as required for MAG_DEPTHS
    // Q: where is maxbin2 noclass? All filtered out already?
    ch_binning_results_gunzipped
            .map { meta, results ->
                def meta_new = meta.clone()
                [ [ 'id': meta_new['id'], 'group': meta_new['group'], 'single_end': meta_new['single_end'], 'assembler': meta_new['assembler'] ], results ]
            }
        .groupTuple (by: 0 )
        .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0 )
        .set { ch_depth_input }

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

    // Group final binned contigs per sample for final output
    ch_binning_results_gunzipped
        .groupTuple(by: 0)
        .set{ ch_binning_results_final }

    emit:
    bins                                         = ch_binning_results_final
    unbinned                                     = SPLIT_FASTA.out.unbinned
    tooshort                                     = METABAT2_METABAT2.out.tooshort
    lowdepth                                     = METABAT2_METABAT2.out.lowdepth
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    metabat2_version                             = METABAT2_METABAT2.out.versions
    metabat2_jgisummarizebamcontigdepths_version = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions
}
