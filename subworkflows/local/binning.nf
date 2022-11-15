/*
 * Binning with MetaBAT2 and MaxBin2
 */

params.mag_depths_options                           = [:]
params.mag_depths_plot_options                      = [:]
params.mag_depths_summary_options                   = [:]

include { METABAT2_METABAT2                     } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                               } from '../../modules/nf-core/maxbin2/main'
include { GUNZIP as GUNZIP_BINS                 } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS               } from '../../modules/nf-core/gunzip/main'

include { CONVERT_DEPTHS                        } from '../../modules/local/convert_depths'
include { ADJUST_MAXBIN2_EXT                    } from '../../modules/local/adjust_maxbin2_ext'
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

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, depths ]
        }

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    // combine depths back with assemblies
    ch_metabat2_input = assemblies
        .map { meta, assembly, bams, bais ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, assembly, bams, bais ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }

    // convert metabat2 depth files to maxbin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( ch_metabat2_input )
        ch_maxbin2_input = CONVERT_DEPTHS.out.output
            .map { meta, assembly, reads, depth ->
                    def meta_new = meta.clone()
                    meta_new['binner'] = 'MaxBin2'

                [ meta_new, assembly, reads, depth ]
            }
        ch_versions = ch_versions.mix(CONVERT_DEPTHS.out.versions.first())
    }

    // main bins for decompressing for MAG_DEPTHS
    ch_final_bins_for_gunzip = Channel.empty()
    // final gzipped bins
    ch_binning_results_gzipped_final = Channel.empty()
    // run binning
    if ( !params.skip_metabat2 ) {
        METABAT2_METABAT2 ( ch_metabat2_input )
        // before decompressing first have to separate and re-group due to limitation of GUNZIP module
        ch_final_bins_for_gunzip = ch_final_bins_for_gunzip.mix( METABAT2_METABAT2.out.fasta.transpose() )
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix( METABAT2_METABAT2.out.fasta )
        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    }
    if ( !params.skip_maxbin2 ) {
        MAXBIN2 ( ch_maxbin2_input )
        ADJUST_MAXBIN2_EXT ( MAXBIN2.out.binned_fastas )
        ch_final_bins_for_gunzip = ch_final_bins_for_gunzip.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins.transpose() )
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins )
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
    // large unbinned contigs from SPLIT_FASTA for decompressing for MAG_DEPTHS,
    // first have to separate and re-group due to limitation of GUNZIP module
    ch_split_fasta_results_transposed = SPLIT_FASTA.out.unbinned.transpose()
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    GUNZIP_BINS ( ch_final_bins_for_gunzip )
    ch_binning_results_gunzipped = GUNZIP_BINS.out.gunzip
    ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions.first())

    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip
    ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions.first())

    // Compute bin depths for different samples (according to `binning_map_mode`)
    // Have to remove binner meta before joining with according depths files,
    // as required for MAG_DEPTHS, but we can add  'binner'
    // info again based on file name and finally group by
    // 'assembler', 'id', 'binner'
    ch_depth_input = ch_binning_results_gunzipped
        .mix(ch_splitfasta_results_gunzipped )
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.remove('binner')
            [ meta_new, bin ]
        }
        .groupTuple (by: 0 )
        .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0 )
        .transpose()
        .map { meta, bin, contig_depths_file ->
            def meta_new = meta.clone()
            meta_new['binner'] = bin.name.split("-")[1]
            [ meta_new, bin, contig_depths_file ]
        }
        .groupTuple (by: [0,2] )

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

    // Group final binned contigs per sample for final output
    ch_binning_results_gunzipped_final = ch_binning_results_gunzipped.groupTuple(by: 0)
    ch_binning_results_gzipped_final   = ch_binning_results_gzipped_final.groupTuple(by: 0)

    emit:
    bins                                         = ch_binning_results_gunzipped_final
    bins_gz                                      = ch_binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped.groupTuple()
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    metabat2depths                               = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    versions                                     = ch_versions
}
