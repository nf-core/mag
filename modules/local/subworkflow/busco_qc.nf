/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

params.busco_db_options      = [:]
params.busco_options         = [:]
params.busco_plot_options    = [:]
params.busco_summary_options = [:]

include { BUSCO_DB_PREPARATION                    } from '../process/busco_db_preparation'        addParams( options: params.busco_db_options      )
include { BUSCO                                   } from '../process/busco'                       addParams( options: params.busco_options         )
include { BUSCO_PLOT                              } from '../process/busco_plot'                  addParams( options: params.busco_plot_options    )
include { BUSCO_SUMMARY                           } from '../process/busco_summary'               addParams( options: params.busco_summary_options )

workflow BUSCO_QC {
    take:
    busco_db_file           // channel: path
    bins                    // channel: [ val(assembler), val(name), path(bin) ]

    main:
    BUSCO_DB_PREPARATION ( busco_db_file )
    BUSCO (
        bins,
        BUSCO_DB_PREPARATION.out.db.collect()
    )

    // group by assembler and sample name for plotting
    BUSCO_PLOT ( BUSCO.out.summary.groupTuple(by: [0,1]) )
    BUSCO_SUMMARY ( BUSCO.out.summary.map{it[2]}.collect() )

    emit:
    summary = BUSCO_SUMMARY.out
    multiqc = BUSCO.out.summary.map{it[2]}
    version = BUSCO.out.version
}
