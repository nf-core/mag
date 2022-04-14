/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */

params.busco_db_options            = [:]
params.busco_options               = [:]
params.busco_save_download_options = [:]
params.busco_plot_options          = [:]
params.busco_summary_options       = [:]

include { BUSCO_DB_PREPARATION            } from '../../modules/local/busco_db_preparation'        addParams( options: params.busco_db_options            )
include { BUSCO                           } from '../../modules/local/busco'                       addParams( options: params.busco_options               )
include { BUSCO_SAVE_DOWNLOAD             } from '../../modules/local/busco_save_download'         addParams( options: params.busco_save_download_options )
include { BUSCO_PLOT                      } from '../../modules/local/busco_plot'                  addParams( options: params.busco_plot_options          )
include { BUSCO_SUMMARY                   } from '../../modules/local/busco_summary'               addParams( options: params.busco_summary_options       )

workflow BUSCO_QC {
    take:
    busco_db_file           // channel: path
    busco_download_folder   // channel: path
    bins                    // channel: [ val(meta), path(bin) ]

    main:
    if (params.busco_reference){
        BUSCO_DB_PREPARATION ( busco_db_file )
        ch_busco_db = BUSCO_DB_PREPARATION.out.db
    } else {
        ch_busco_db = Channel.empty()
    }
    BUSCO (
        bins,
        ch_busco_db.collect().ifEmpty([]),
        busco_download_folder.collect().ifEmpty([])
    )
    if (params.save_busco_reference){
        // publish files downloaded by Busco
        ch_downloads = BUSCO.out.busco_downloads.groupTuple().map{lin,downloads -> downloads[0]}.toSortedList().flatten()
        BUSCO_SAVE_DOWNLOAD ( ch_downloads )
    }
    // group by assembler, binner and sample name for plotting
    // note: failed bins and bins with no domain, i.e. with viral lineages selected, will not be represented in these plots
    ch_results_busco_plot = BUSCO.out.summary_specific.groupTuple(by: 0)
    if (!params.busco_reference){
        ch_results_busco_plot = ch_results_busco_plot.mix(BUSCO.out.summary_domain.groupTuple(by: 0))
    }
    BUSCO_PLOT ( ch_results_busco_plot )

    BUSCO_SUMMARY (
        BUSCO.out.summary_domain.map{it[1]}.collect().ifEmpty([]),
        BUSCO.out.summary_specific.map{it[1]}.collect().ifEmpty([]),
        BUSCO.out.failed_bin.map{it[1]}.collect().ifEmpty([])
    )

    emit:
    summary     = BUSCO_SUMMARY.out.summary
    failed_bin  = BUSCO.out.failed_bin.map{it[1]}
    multiqc     = BUSCO.out.summary_domain.map{it[1]}
    versions    = BUSCO.out.versions
}
