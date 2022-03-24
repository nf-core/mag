/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { DASTOOL_SCAFFOLDS2BIN as DASTOOL_SCAFFOLDS2BIN_METABAT2 } from '../../modules/nf-core/modules/dastool/scaffolds2bin/main.nf'
include { DASTOOL_SCAFFOLDS2BIN as DASTOOL_SCAFFOLDS2BIN_MAXBIN2  } from '../../modules/nf-core/modules/dastool/scaffolds2bin/main.nf'
include { DASTOOL_DASTOOL                                         } from '../../modules/nf-core/modules/dastool/dastool/main.nf'

workflow BINNING_REFINEMENT {
    take:
    contigs        //
    bins           // channel: [ val(meta), path(bins) ]

    main:
    ch_versions = Channel.empty()

    ch_contigs_for_dastool = contigs
                                .map {
                                    meta, assembly, bams, bais ->
                                        def meta_new = meta.clone()
                                        [ meta_new, assembly ]
                                }
                                .dump(tag: "binrefine_contigs")
    bins.dump(tag: "binrefine_bins")

    ch_bins_for_scaffolds2bin = bins
                                    .branch {
                                        metabat2: it[0]['binner'] == 'MetaBAT2'
                                        maxbin2:  it[0]['binner'] == 'MaxBin2'
                                    }

    // Run on each bin file separately
    DASTOOL_SCAFFOLDS2BIN_METABAT2 ( ch_bins_for_scaffolds2bin.metabat2, "fa")
    DASTOOL_SCAFFOLDS2BIN_MAXBIN2  ( ch_bins_for_scaffolds2bin.maxbin2, "fasta")

    ch_versions = ch_versions.mix(DASTOOL_SCAFFOLDS2BIN_METABAT2.out.versions)
    ch_versions = ch_versions.mix(DASTOOL_SCAFFOLDS2BIN_MAXBIN2.out.versions)

    // Concatenate each binner table per sample (based on contig names in FASTA file, not fasta file itself, so should be OK)

    // Run DAStool
    //DASTOOL_DASTOOL()
    //ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions)

    emit:
    versions               = ch_versions
}
