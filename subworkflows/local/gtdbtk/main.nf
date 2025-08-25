/*
 * GTDB-Tk bin classification, using BUSCO QC to filter bins
 */

include { GTDBTK_CLASSIFYWF     } from '../../../modules/nf-core/gtdbtk/classifywf/main'

include { GTDBTK_DB_PREPARATION } from '../../../modules/local/gtdbtk_db_preparation/main'
include { GTDBTK_SUMMARY        } from '../../../modules/local/gtdbtk_summary/main'

workflow GTDBTK {
    take:
    ch_bins           // channel: [ val(meta), [bins] ]
    ch_bin_qc_summary // channel: path
    val_gtdb          // value: path
    val_gtdb_mash     // value: path

    main:
    ch_versions = Channel.empty()

    // Collect bin quality metrics
    qc_columns = [
        busco: ['Input_file', 'Complete', 'Duplicated'],
        checkm: ['Bin Id', 'Completeness', 'Contamination'],
        checkm2: ['Name', 'Completeness', 'Contamination'],
    ]

    ch_bin_metrics = ch_bin_qc_summary
        .splitCsv(header: true, sep: '\t')
        .map { row -> qc_columns[params.binqc_tool].collect { col -> row[col] } }
        .filter { row -> row[1] != '' }
        .map { row ->
            row = [row[0]] + row[1..2].collect { value -> "${value}".toDouble() }
            // CheckM / CheckM2 removes the .fa extension from the bin name
            if (params.binqc_tool in ['checkm', 'checkm2']) {
                row[0] = row[0] + '.fa'
            }
            row
        }

    // Filter bins based on collected metrics: completeness, contamination
    ch_filtered_bins = ch_bins
        .transpose()
        .map { meta, bin -> [bin.getName(), bin, meta] }
        .join(ch_bin_metrics, failOnDuplicate: true)
        .map { _bin_name, bin, meta, completeness, contamination -> [meta, bin, completeness, contamination] }
        .branch {
            passed: (it[2] != -1 && it[2] >= params.gtdbtk_min_completeness && it[3] != -1 && it[3] <= params.gtdbtk_max_contamination)
            return [it[0], it[1]]
            discarded: (it[2] == -1 || it[2] < params.gtdbtk_min_completeness || it[3] == -1 || it[3] > params.gtdbtk_max_contamination)
            return [it[0], it[1]]
        }

    if (val_gtdb.extension == 'gz') {
        // Expects to be tar.gz!
        ch_db_for_gtdbtk = GTDBTK_DB_PREPARATION(val_gtdb).db
        ch_versions = ch_versions.mix(GTDBTK_DB_PREPARATION.out.versions)
    }
    else if (val_gtdb.isDirectory()) {
        ch_db_for_gtdbtk = [val_gtdb.simpleName, val_gtdb]
    }
    else {
        error("Unsupported object given to --gtdb, database must be supplied as either a directory or a .tar.gz file!")
    }


    // Print warning why GTDB-TK summary empty if passed channel gets no files
    ch_filtered_bins.passed
        .count()
        .map { it == 0 ? log.warn("No contigs passed GTDB-TK min. completeness filters. GTDB-TK summary will execute but results will be empty!") : "" }


    GTDBTK_CLASSIFYWF(
        ch_filtered_bins.passed.groupTuple(),
        ch_db_for_gtdbtk,
        params.gtdbtk_pplacer_useram ? false : true,
        val_gtdb_mash,
    )
    ch_versions = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions)

    GTDBTK_SUMMARY(
        ch_filtered_bins.discarded.map { it[1] }.collect().ifEmpty([]),
        GTDBTK_CLASSIFYWF.out.summary.map { it[1] }.collect().ifEmpty([]),
        [],
        [],
    )
    ch_versions = ch_versions.mix(GTDBTK_SUMMARY.out.versions)

    emit:
    summary       = GTDBTK_SUMMARY.out.summary
    multiqc_files = GTDBTK_SUMMARY.out.summary
    versions      = ch_versions
}
