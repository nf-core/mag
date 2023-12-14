/*
 * GTDB-Tk bin classification, using BUSCO QC to filter bins
 */

include { GTDBTK_DB_PREPARATION } from '../../modules/local/gtdbtk_db_preparation'
include { GTDBTK_CLASSIFYWF     } from '../../modules/nf-core/gtdbtk/classifywf/main'
include { GTDBTK_SUMMARY        } from '../../modules/local/gtdbtk_summary'

workflow GTDBTK {
    take:
    bins              // channel: [ val(meta), [bins] ]
    busco_summary     // channel: path
    checkm_summary    // channel: path
    gtdb              // channel: path
    gtdb_mash         // channel: path

    main:
    // Filter bins: classify only medium & high quality MAGs
    ch_bin_metrics = Channel.empty()
    if ( params.binqc_tool == 'busco' ){
        // Collect completeness and contamination metrics from busco summary
        ch_bin_metrics = busco_summary
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                        def completeness  = -1
                        def contamination = -1
                        def missing, duplicated
                        if (params.busco_db.getBaseName().contains('odb10')) {
                            missing    = row.'%Missing (specific)'      // TODO or just take '%Complete'?
                            duplicated = row.'%Complete and duplicated (specific)'
                        } else {
                            missing    = row.'%Missing (domain)'
                            duplicated = row.'%Complete and duplicated (domain)'
                        }
                        if (missing != '') completeness = 100.0 - Double.parseDouble(missing)
                        if (duplicated != '') contamination = Double.parseDouble(duplicated)
                        [row.'GenomeBin', completeness, contamination]
            }
    } else {
        // Collect completeness and contamination metrics from checkm summary
        ch_bin_metrics = checkm_summary
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                        def completeness  = Double.parseDouble(row.'Completeness')
                        def contamination = Double.parseDouble(row.'Contamination')
                        [row.'Bin Id' + ".fa", completeness, contamination]
            }
    }


    // Filter bins based on collected metrics: completeness, contamination
    ch_filtered_bins = bins
        .transpose()
        .map { meta, bin -> [bin.getName(), bin, meta]}
        .join(ch_bin_metrics, failOnDuplicate: true, failOnMismatch: true)
        .map { bin_name, bin, meta, completeness, contamination -> [meta, bin, completeness, contamination] }
        .branch {
            passed: (it[2] != -1 && it[2] >= params.gtdbtk_min_completeness && it[3] != -1 && it[3] <= params.gtdbtk_max_contamination)
                return [it[0], it[1]]
            discarded: (it[2] == -1 || it[2] < params.gtdbtk_min_completeness || it[3] == -1 || it[3] > params.gtdbtk_max_contamination)
                return [it[0], it[1]]
        }

    if ( gtdb.extension == 'gz' ) {
        // Expects to be tar.gz!
        ch_db_for_gtdbtk = GTDBTK_DB_PREPARATION ( gtdb ).db
    } else if ( gtdb.isDirectory() ) {
        // Make up meta id to match expected channel cardinality for GTDBTK
        ch_db_for_gtdbtk = Channel
                            .of(gtdb)
                            .map{
                                [ it.toString().split('/').last(), it ]
                            }
                            .collect()
    } else {
        error("Unsupported object given to --gtdb, database must be supplied as either a directory or a .tar.gz file!")
    }

    GTDBTK_CLASSIFYWF (
        ch_filtered_bins.passed.groupTuple(),
        ch_db_for_gtdbtk,
        gtdb_mash
    )

    GTDBTK_SUMMARY (
        ch_filtered_bins.discarded.map{it[1]}.collect().ifEmpty([]),
        GTDBTK_CLASSIFYWF.out.summary.collect().ifEmpty([]),
        [],
        // GTDBTK_CLASSIFYWF.out.filtered.collect().ifEmpty([]),
        []
        // GTDBTK_CLASSIFYWF.out.failed.collect().ifEmpty([])
    )

    emit:
    gtdbtk_summaries    = GTDBTK_CLASSIFYWF.out.summary
    summary             = GTDBTK_SUMMARY.out.summary
    versions            = GTDBTK_CLASSIFYWF.out.versions
}
