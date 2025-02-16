/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { METABAT2_METABAT2                                            } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS                         } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                                                      } from '../../modules/nf-core/maxbin2/main'
include { GUNZIP as GUNZIP_BINS                                        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS                                      } from '../../modules/nf-core/gunzip/main'

include { CONVERT_DEPTHS                        } from '../../modules/local/convert_depths'
include { ADJUST_MAXBIN2_EXT                    } from '../../modules/local/adjust_maxbin2_ext'
include { SPLIT_FASTA                           } from '../../modules/local/split_fasta'
include { FASTA_BINNING_CONCOCT                 } from '../../subworkflows/nf-core/fasta_binning_concoct/main'

workflow BINNING {
    take:
    assemblies           // channel: [ val(meta), path(assembly), path(bams), path(bais) ]

    main:

    ch_versions = Channel.empty()

    // generate coverage depths for each contig
    ch_summarizedepth_input = assemblies
                                .map { meta, _assembly, bams, bais ->
                                    [ meta, bams, bais ]
                                }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, depths ]
        }

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    // combine depths back with assemblies
    ch_metabat2_input = assemblies
        .map { meta, assembly, bams, bais ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, assembly, bams, bais ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, assembly, _bams, _bais, depths ->
            [ meta, assembly, depths ]
        }

    // convert metabat2 depth files to maxbin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( ch_metabat2_input )
        ch_maxbin2_input = CONVERT_DEPTHS.out.output
            .map { meta, assembly, reads, depth ->
                    def meta_new = meta + [binner: 'MaxBin2']
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
    if ( !params.skip_concoct ){

        ch_concoct_input = assemblies
                            .map { meta, bins, bams, bais ->
                                def meta_new = meta + [binner: 'CONCOCT']
                                [ meta_new, bins, bams, bais ]
                            }
                            .multiMap {
                                meta, bins, bams, bais ->
                                    bins: [ meta, bins ]
                                    bams: [ meta, bams, bais ]
                            }

        FASTA_BINNING_CONCOCT ( ch_concoct_input.bins, ch_concoct_input.bams )
        ch_final_bins_for_gunzip = ch_final_bins_for_gunzip.mix( FASTA_BINNING_CONCOCT.out.bins.transpose() )
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix( FASTA_BINNING_CONCOCT.out.bins )
        ch_versions = ch_versions.mix(FASTA_BINNING_CONCOCT.out.versions)
    }

    // decide which unbinned fasta files to further filter, depending on which binners selected
    // NOTE: CONCOCT does not produce 'unbins' itself, therefore not included here.
    if ( !params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = METABAT2_METABAT2.out.unbinned
    } else if ( params.skip_metabat2 && !params.skip_maxbin2 ) {
        ch_input_splitfasta = MAXBIN2.out.unbinned_fasta
    } else if ( params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = Channel.empty()
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
        .groupTuple(by: 0)

    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip
        .groupTuple(by: 0)

    ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions.first())
    ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions.first())

    emit:
    bins                                         = ch_binning_results_gunzipped
    bins_gz                                      = ch_binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    metabat2depths                               = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    versions                                     = ch_versions
}
