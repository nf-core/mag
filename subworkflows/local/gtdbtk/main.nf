/*
 * GTDB-Tk bin classification, using BUSCO QC to filter bins
 */

include { GTDBTK_CLASSIFYWF     } from '../../../modules/nf-core/gtdbtk/classifywf/main'

include { GTDBTK_DB_PREPARATION } from '../../../modules/local/gtdbtk_db_preparation/main'
include { GTDBTK_SUMMARY        } from '../../../modules/local/gtdbtk_summary/main'

workflow GTDBTK {
    take:
    ch_bins           // [val(meta), path(fasta)]
    ch_bin_qc_summary // path
    val_gtdb          // path

    main:
    ch_versions = channel.empty()

    // Collect bin quality metrics
    qc_columns = [
        busco: ['bin_qc_tool', 'Input_file', 'Complete', 'Duplicated'],
        checkm: ['bin_qc_tool', 'Bin Id', 'Completeness', 'Contamination'],
        checkm2: ['bin_qc_tool', 'Name', 'Completeness', 'Contamination'],
    ]

    // 1. Take the pre-split CSV, and select just the relevant columns
    // 2. Reorder columns to: bin name, tool, completeness, contamination/duplication
    // 3. Convert completeness/contamination to double for consistency
    //4. add .fa back to bin name if missing
    ch_bin_metrics = ch_bin_qc_summary
        .map { row -> qc_columns[row.bin_qc_tool].collect { col -> row[col] } }
        .map { row ->
            // Initial order
            // row[0] = bin_qc_tool
            // row[1] = bin name
            // row[2] = completeness
            // row[3] = contamination (Checkm*)/duplication (busco)
            def row_reordered = [row[1], row[0]] + row[2..3].collect { value -> "${value}".toDouble() }
            // CheckM / CheckM2 removes the .fa extension from the bin name
            if (row_reordered[1] in ['checkm', 'checkm2']) {
                row_reordered[0] = row_reordered[0] + '.fa'
            }
            return row_reordered
        }

    // Filter bins based on collected metrics: completeness, contamination
    // 1. Generate key name from bin file name
    // 2. Join with metrics (but drop bin_name)
    // 3. Group all metrics files per bin (in case multiple QC tools were used)
    // 3. Branch based on completeness/contamination
    ch_filtered_bins = ch_bins
        .transpose()
        .map { meta, bin -> [bin.getName(), bin, meta] }
        .join(ch_bin_metrics)
        .map { bin_name, bin, meta, _bin_qc_tool, completeness, contamination ->
            [bin_name, meta, bin, completeness, contamination]
        }
        .groupTuple(by: 0)
        .branch { _bin_name, meta, bin, completeness, contamination ->
            passed: (
                completeness.any { bin_completeness -> bin_completeness != -1 } &&
                completeness.any { bin_completeness -> bin_completeness >= params.gtdbtk_min_completeness } &&
                contamination.any { bin_contamination -> bin_contamination != -1 } &&
                contamination.any { bin_contamination -> bin_contamination <= params.gtdbtk_max_contamination }
            )
                return [meta[0], bin[0]]
            discarded: true
                return [meta[0], bin[0]]
        }
    // Note we have to call `meta[0], bin[0]` because of the groupTuple above

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

    GTDBTK_CLASSIFYWF(
        ch_filtered_bins.passed.groupTuple(),
        ch_db_for_gtdbtk,
        params.gtdbtk_pplacer_useram ? false : true,
    )
    ch_versions = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions)

    // Print warning why GTDB-TK summary empty if passed channel gets no files
    ch_filtered_bins.passed
        .count()
        .combine(ch_filtered_bins.discarded.count())
        .subscribe { passed, failed ->
            if ((passed + failed) > 0 && passed == 0) {
                log.warn("No contigs passed GTDB-TK min. completeness filters. GTDB-TK summary will execute but results will be empty!")
            }
        }

    GTDBTK_SUMMARY(
        ch_filtered_bins.discarded.map { _meta, bin -> bin }.collect().ifEmpty([]),
        GTDBTK_CLASSIFYWF.out.summary.map { _meta, summary -> summary }.collect().ifEmpty { ([]) },
        [],
        [],
    )
    ch_versions = ch_versions.mix(GTDBTK_SUMMARY.out.versions)

    emit:
    summary       = GTDBTK_SUMMARY.out.summary
    multiqc_files = GTDBTK_SUMMARY.out.summary
    versions      = ch_versions
}
