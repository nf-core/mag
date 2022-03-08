/*
 * Bin refinement with DASTool
 */

include { DASTOOL_SCAFFOLDS2BIN as DASTOOL_SCAFFOLDS2BIN_METABAT2  } from '../../modules/nf-core/modules/dastool/scaffolds2bin/main.nf'
include { DASTOOL_SCAFFOLDS2BIN as DASTOOL_SCAFFOLDS2BIN_MAXBIN2   } from '../../modules/nf-core/modules/dastool/scaffolds2bin/main.nf'
include { DASTOOL_DASTOOL                                          } from '../../modules/nf-core/modules/dastool/dastool/main.nf'

workflow BINNING_REFINEMENT {
    take:
    contigs        // channel: [ val(meta), path(fasta)) ]
    bins           // channel: [ val(meta), path(binned_fasta_gz) ]

    main:
    ch_versions = Channel.empty()

    // run on each bin file separately due to different suffixes
    bins
        .dump(tag: "input_to_binrefinement")
        .branch{
            maxbin2: it[0]['binner'] == 'MaxBin2'
            metabat2: it[0]['binner'] == 'MetaBAT2'
        }
        .set{ ch_for_scaffolds2bin }

    // generate bin files
    if ( !params.skip_metabat2 & !params.skip_maxbin2 ) {
        DASTOOL_SCAFFOLDS2BIN_METABAT2 ( ch_for_scaffolds2bin.metabat2, "fa")
        DASTOOL_SCAFFOLDS2BIN_MAXBIN2 ( ch_for_scaffolds2bin.maxbin2, "fasta")
        DASTOOL_SCAFFOLDS2BIN_METABAT2.out.scaffolds2bin
            .mix(DASTOOL_SCAFFOLDS2BIN_MAXBIN2.out.scaffolds2bin)
            .set { ch_s2b_for_dastool }
        ch_versions = ch_versions.mix( DASTOOL_SCAFFOLDS2BIN_METABAT2.out.versions, DASTOOL_SCAFFOLDS2BIN_MAXBIN2.out.versions )
    } else if ( !params.skip_metabat2 & params.skip_maxbin2 ){
        DASTOOL_SCAFFOLDS2BIN_METABAT2 ( ch_for_scaffolds2bin.metabat2, "fa")
        DASTOOL_SCAFFOLDS2BIN_METABAT2.out.scaffolds2bin.set { ch_s2b_for_dastool }
        ch_versions = ch_versions.mix( DASTOOL_SCAFFOLDS2BIN_METABAT2.out.versions )
    } else if ( params.skip_metabat2 & !params.skip_maxbin2 ){
        DASTOOL_SCAFFOLDS2BIN_MAXBIN2 ( ch_for_scaffolds2bin.maxbin2, "fasta")
        DASTOOL_SCAFFOLDS2BIN_MAXBIN2.out.scaffolds2bin.set { ch_s2b_for_dastool }
        ch_versions = ch_versions.mix( DASTOOL_SCAFFOLDS2BIN_MAXBIN2.out.versions )
    }

    // combine per sample + assembler combination
    // TODO scaffolds2bin remove binner from meta, group
    // TODO contigs remove binner from meta, group
    // mergin into input for DASTool
    ch_s2b_for_dastool
        .dump(tag: "s2b_prior_meta_reduc")
        .map { meta, results ->
            def meta_new = meta.clone()
            [ [ 'id': meta_new['id'], 'group': meta_new['group'], 'single_end': meta_new['single_end'], 'assembler': meta_new['assembler'] ], results ]
        }
        .dump(tag: "s2b_prior_merge")
        .groupTuple()
        .set{ ch_s2b_for_join }

    contigs.dump(tag: "input_contigs_to_s2b_join").join(ch_s2b_for_join, by: 0 ).set{ ch_input_for_dastool }

    // TODO now picking up noclass/unbinned?
    ch_input_for_dastool.dump(tag: "input_to_dastool")

    // Run DAS_Tool
    // TODO offer extra parameters
    DASTOOL_DASTOOL( ch_input_for_dastool, [], [], [] )


    // Concatenate each binner table per sample (based on contig names in FASTA file, not fasta file itself, so should be OK)



    emit:
    //refined_bins
    versions = ch_versions

}
