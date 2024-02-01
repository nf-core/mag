/*
 * GUNC: Detection and quantification of genome chimerism based on lineage homogeneity
 */

include { GUNC_DOWNLOADDB   } from '../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN          } from '../../modules/nf-core/gunc/run/main'
include { GUNC_MERGECHECKM  } from '../../modules/nf-core/gunc/mergecheckm/main'

workflow GUNC_QC {
    take:
    ch_bins         // [ [ meta] , fasta           ], input bins (mandatory)
    ch_gunc_db      // [ db                        ], presupplied GUNC database (optional)
    ch_checkm_table // [ [ meta ], checkm_qa_table ], extended checkm table from CHECKM_QA, (optional)

    main:
    ch_versions = Channel.empty()

    if ( params.gunc_db ) {
        ch_db_for_gunc = ch_gunc_db
    } else {
        ch_db_for_gunc = GUNC_DOWNLOADDB( params.gunc_database_type ).db
        ch_versions.mix( GUNC_DOWNLOADDB.out.versions )
    }


    GUNC_RUN ( ch_bins, ch_db_for_gunc )
    ch_versions.mix( GUNC_RUN.out.versions )

    // Make sure to keep directory in sync with modules.conf
    GUNC_RUN.out.maxcss_level_tsv
        .map{it[1]}
        .collectFile(name: "gunc_summary.tsv", keepHeader: true, storeDir: "${params.outdir}/GenomeBinning/QC/")

    if ( params.binqc_tool == 'checkm' ) {

        ch_input_to_mergecheckm = GUNC_RUN.out.maxcss_level_tsv
                                    .combine(ch_checkm_table, by: 0)

        GUNC_MERGECHECKM ( ch_input_to_mergecheckm )
        ch_versions.mix( GUNC_MERGECHECKM.out.versions )

    // Make sure to keep directory in sync with modules.conf
    GUNC_MERGECHECKM.out.tsv
        .map{it[1]}
        .collectFile(name: "gunc_checkm_summary.tsv", keepHeader: true, storeDir: "${params.outdir}/GenomeBinning/QC/")
    }

    emit:
    versions = ch_versions

}
