/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

params.busco_db_options      = [:]
params.busco_options         = [:]
params.busco_plot_options    = [:]
params.busco_summary_options = [:]

include { BUSCO_DB_PREPARATION                    } from '../../modules/local/busco_db_preparation'        addParams( options: params.busco_db_options      )
include { BUSCO                                   } from '../../modules/local/busco'                       addParams( options: params.busco_options         )
include { BUSCO_PLOT                              } from '../../modules/local/busco_plot'                  addParams( options: params.busco_plot_options    )
include { BUSCO_SUMMARY                           } from '../../modules/local/busco_summary'               addParams( options: params.busco_summary_options )

workflow BUSCO_QC {
    take:
    busco_db_file           // channel: path
    bins                    // channel: [ val(meta), path(bin) ]

    main:
    BUSCO_DB_PREPARATION ( busco_db_file )
    BUSCO (
        bins,
        BUSCO_DB_PREPARATION.out.db
    )

    // group by assembler and sample name for plotting
    BUSCO_PLOT ( BUSCO.out.summary.groupTuple(by: 0) )
    BUSCO_SUMMARY ( BUSCO.out.summary.map{it[1]}.collect() )

    emit:
    summary = BUSCO_SUMMARY.out
    multiqc = BUSCO.out.summary.map{it[1]}
    version = BUSCO.out.version
}
