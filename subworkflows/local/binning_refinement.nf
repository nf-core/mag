/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { DASTOOL_SCAFFOLDS2BIN } from '../../modules/dastool/scaffolds2bin/main.nf'
include { DASTOOL_DASTOOL       } from '../../modules/dastool/dastool/main.nf'

workflow BINNING {
    take:
    contigs        //
    bins           // channel: [ val(meta), path(bins) ]

    main:


    // Run on each bin file separately
    DASTOOL_SCAFFOLDS2BIN ( METABAT2_METABAT2.out.fasta.collect(), "fa")

    // Concatenate each binner table per sample (based on contig names in FASTA file, not fasta file itself, so should be OK)

    // Run DAStool



    emit:
    bins                                         = ch_binning_results_final
    unbinned                                     = SPLIT_FASTA.out.unbinned
    tooshort                                     = METABAT2_METABAT2.out.tooshort
    lowdepth                                     = METABAT2_METABAT2.out.lowdepth
    depths_summary                               = MAG_DEPTHS_SUMMARY.out.summary
    metabat2_version                             = METABAT2_METABAT2.out.versions
    metabat2_jgisummarizebamcontigdepths_version = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions
    maxbin2_version                              = MAXBIN2.out.versions
}
