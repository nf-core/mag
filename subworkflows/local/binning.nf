/*
 * Binning with MetaBAT2 and MaxBin2
 */

params.mag_depths_options                           = [:]
params.mag_depths_plot_options                      = [:]
params.mag_depths_summary_options                   = [:]

include { METABAT2_METABAT2                     } from '../../modules/nf-core/modules/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../modules/nf-core/modules/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                               } from '../../modules/nf-core/modules/maxbin2/main'
include { GUNZIP as GUNZIP_BINS                 } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS               } from '../../modules/nf-core/modules/gunzip/main'

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

    ch_versions = Channel.empty()

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

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions)

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
        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions)
    }

    if ( !params.skip_maxbin2 ) {
        MAXBIN2 ( ch_maxbin2_input )
        ch_versions = ch_versions.mix(MAXBIN2.out.versions)
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
    // MAG_DEPTHS, first have to separate and re-group due to limitation of
    // GUNZIP module
    METABAT2_METABAT2.out.fasta.transpose().set { ch_metabat2_results_transposed }
    MAXBIN2.out.binned_fastas.transpose().set { ch_maxbin2_results_transposed }

    SPLIT_FASTA.out.unbinned.transpose().set { ch_split_fasta_results_transposed }
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    // mix all bins together
    ch_metabat2_results_transposed
        .mix( ch_maxbin2_results_transposed )
        .set { ch_final_bins_for_gunzip }

    GUNZIP_BINS ( ch_final_bins_for_gunzip )
    GUNZIP_BINS.out.gunzip
        .set{ ch_binning_results_gunzipped }
    ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions)

    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
    GUNZIP_UNBINS.out.gunzip
        .set{ ch_splitfasta_results_gunzipped }
    ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions)

    // Compute bin depths for different samples (according to `binning_map_mode`)
    // Have to remove binner meta for grouping to mix back with original depth
    // files, as required for MAG_DEPTHS
    ch_binning_results_gunzipped
        .mix(ch_splitfasta_results_gunzipped )
        .map { meta, results ->
            def meta_new = meta.clone()
            [ [ 'id': meta_new['id'], 'group': meta_new['group'], 'single_end': meta_new['single_end'], 'assembler': meta_new['assembler'] ], results ]
        }
        .groupTuple (by: 0 )
        .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0 )
        .set { ch_depth_input }

    MAG_DEPTHS ( ch_depth_input )
    ch_versions = ch_versions.mix(MAG_DEPTHS.out.versions)

    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    // create file containing group information for all samples
    ch_sample_groups = reads
        .collectFile(name:'sample_groups.tsv'){ meta, reads -> meta.id + '\t' + meta.group + '\n' }

    // filter MAG depth files: use only those for plotting that contain depths for > 2 samples
    ch_mag_depths_plot = MAG_DEPTHS.out.depths
        .map { meta, depth_file -> if (getColNo(depth_file) > 2) [meta, depth_file] }

    MAG_DEPTHS_PLOT ( ch_mag_depths_plot, ch_sample_groups.collect() )
    MAG_DEPTHS_SUMMARY ( MAG_DEPTHS.out.depths.map{it[1]}.collect() )
    ch_versions = ch_versions.mix(MAG_DEPTHS_PLOT.out.versions)
    ch_versions = ch_versions.mix(MAG_DEPTHS_SUMMARY.out.versions)

    // Group final binned contigs per sample for final output
    ch_binning_results_gunzipped
        .groupTuple(by: 0)
        .set{ ch_binning_results_gunzipped_final }

    METABAT2_METABAT2.out.fasta.mix(MAXBIN2.out.binned_fastas)
        .groupTuple(by: 0)
        .set{ ch_binning_results_gzipped_final }

    SPLIT_FASTA.out.unbinned

    emit:
    bins                                         = ch_binning_results_gunzipped_final
    bins_gz                                      = ch_binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped.groupTuple()
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    versions                                     = ch_versions
}
