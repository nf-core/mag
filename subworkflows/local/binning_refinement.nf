/*
 * Binning with MetaBAT2 and MaxBin2
 */

include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_METABAT2 } from '../../modules/nf-core/modules/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_MAXBIN2  } from '../../modules/nf-core/modules/dastool/fastatocontig2bin/main.nf'
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

    ch_bins_for_fastatocontig2bin = bins
                                    .branch {
                                        metabat2: it[0]['binner'] == 'MetaBAT2'
                                        maxbin2:  it[0]['binner'] == 'MaxBin2'
                                    }

    // Run on each bin file separately
    DASTOOL_FASTATOCONTIG2BIN_METABAT2 ( ch_bins_for_fastatocontig2bin.metabat2, "fa")
    DASTOOL_FASTATOCONTIG2BIN_MAXBIN2 ( ch_bins_for_fastatocontig2bin.maxbin2, "fasta")

    ch_fastatocontig2bin_for_dastool = Channel.empty()
    ch_fastatocontig2bin_for_dastool = ch_fastatocontig2bin_for_dastool
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_METABAT2.out.fastatocontig2bin)
                                    .mix(DASTOOL_FASTATOCONTIG2BIN_MAXBIN2.out.fastatocontig2bin)
                                    .map {
                                        meta, fastatocontig2bin ->
                                            def meta_new = meta.clone()
                                            meta_new.remove('binner')
                                            [ meta_new, fastatocontig2bin ]
                                    }
                                    .groupTuple(by: 0)

    ch_input_for_dastool = ch_contigs_for_dastool.join(ch_fastatocontig2bin_for_dastool, by: 0).dump(tag: "input_to_dastool")

    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_METABAT2.out.versions.first())
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_MAXBIN2.out.versions.first())

    // Concatenate each binner table per sample (based on contig names in FASTA file, not fasta file itself, so should be OK)

    // Run DAStool
    DASTOOL_DASTOOL(ch_input_for_dastool, [], [])
    ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions.first())

    emit:
    dastool     = DASTOOL_DASTOOL.out.bins
    versions    = ch_versions
}
