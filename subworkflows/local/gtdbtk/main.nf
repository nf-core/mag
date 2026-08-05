/*
 * GTDB-Tk bin classification, using BUSCO QC to filter bins
 */

include { GTDBTK_CLASSIFYWF     } from '../../../modules/nf-core/gtdbtk/classifywf/main'

include { GTDBTK_DB_PREPARATION } from '../../../modules/local/gtdbtk_db_preparation/main'
include { GTDBTK_SUMMARY        } from '../../../modules/local/gtdbtk_summary/main'

workflow GTDBTK {
    take:
    ch_bins // [val(meta), [path(fasta)]]  bins grouped per binner-sample
    ch_qc_metrics // [val(meta), val(tool), path(summary)]  per-group, per-tool QC summary
    val_gtdb // path

    main:
    ch_versions = channel.empty()

    // QC summary columns per tool: [bin name, completeness, contamination/duplication]
    qc_columns = [
        busco: ['Input_file', 'Complete', 'Duplicated'],
        checkm: ['Bin Id', 'Completeness', 'Contamination'],
        checkm2: ['Name', 'Completeness', 'Contamination'],
    ]

    // number of QC tools that emit one summary per bin group, used to release each
    // group's metrics as soon as its own QC finishes instead of waiting for all groups
    def n_qc_tools = [params.run_busco, params.run_checkm, params.run_checkm2].count { enabled -> enabled }

    // gather the per-tool summaries for each bin group (non-blocking via groupKey)
    // remainder: true flushes groups whose QC tool failed (fewer summaries than
    // n_qc_tools) at channel close, so a surviving tool can still classify them
    ch_group_metrics = ch_qc_metrics
        .map { meta, tool, summary -> [groupKey(meta, n_qc_tools), meta, tool, summary] }
        .groupTuple(by: 0, remainder: true)
        .map { _key, metas, tools, summaries -> [metas[0], [tools, summaries].transpose()] }

    // Filter bins within each group by completeness/contamination across QC tools.
    ch_filtered_bins = ch_bins
        .join(ch_group_metrics)
        .map { meta, bins, tool_summaries ->
            // build per-bin metrics: bin name -> [ [completeness, contamination], ... ] (one reading per QC tool)
            def metrics = [:]
            tool_summaries.each { tool, summary ->
                def cols = qc_columns[tool]
                summary
                    .splitCsv(header: true, sep: '\t')
                    .each { row ->
                        // CheckM / CheckM2 strip the .fa extension from the bin name
                        def bin_name = tool == 'busco' ? row[cols[0]] : row[cols[0]] + '.fa'
                        def completeness = "${row[cols[1]]}".toDouble()
                        def contamination = "${row[cols[2]]}".toDouble()

                        def readings = metrics.get(bin_name, [])
                        // a negative value means the tool could not assess the bin, so we drop the whole reading
                        if (completeness >= 0 && contamination >= 0) {
                            readings << [completeness: completeness, contamination: contamination]
                        }
                        metrics[bin_name] = readings
                    }
            }

            // drop bins with no QC metric, then split the rest:
            // a bin passes if any single tool clears both thresholds together
            def (passed, discarded) = bins
                .findAll { bin -> metrics[bin.getName() - ~/\.gz$/] != null }
                .split { bin ->
                    metrics[bin.getName() - ~/\.gz$/].any { reading ->
                        reading.completeness >= params.gtdbtk_min_completeness && reading.contamination <= params.gtdbtk_max_contamination
                    }
                }

            [meta, passed, discarded]
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
        error("Unsupported object given to --gtdb_db, database must be supplied as either a directory or a .tar.gz file!")
    }

    // Warn if no QC data was available (all QC tools likely failed)
    ch_qc_metrics
        .count()
        .filter { count -> count == 0 }
        .subscribe { _count ->
            log.warn("[nf-core/mag] No bin QC results were available. Skipping GTDB-Tk classification.")
        }

    ch_passed_bins = ch_filtered_bins
        .map { meta, passed, _discarded -> [meta, passed.toSorted { bin -> bin.name }] }
        .filter { _meta, passed -> passed }
    ch_passed_flat = ch_passed_bins.flatMap { _meta, passed -> passed }
    ch_discarded_bins = ch_filtered_bins.flatMap { _meta, _passed, discarded -> discarded }

    // Optionally combine all bins into a single GTDB-Tk job, instead of one job
    // per bin group, which can be more efficient for very large bin counts
    if (params.gtdbtk_single_job) {
        ch_bins_for_classifywf = ch_passed_flat
            .map { bin ->
                [
                    [id: 'all_bins', assembler: 'all', binner: 'all', domain: 'all', refinement: 'all'],
                    bin,
                ]
            }
            .groupTuple()
            // Sort bins by name so GTDB-Tk receives a deterministic input order.
            .map { meta, bins -> [meta, bins.toSorted { bin -> bin.name }] }
    }
    else {
        ch_bins_for_classifywf = ch_passed_bins
    }

    GTDBTK_CLASSIFYWF(
        ch_bins_for_classifywf,
        ch_db_for_gtdbtk,
        params.gtdbtk_pplacer_useram ? false : true,
    )

    // Print warning why GTDB-TK summary empty if passed channel gets no files
    ch_passed_flat
        .count()
        .combine(ch_discarded_bins.count())
        .subscribe { passed, failed ->
            if ((passed + failed) > 0 && passed == 0) {
                log.warn("[nf-core/mag] No contigs passed GTDB-TK min. completeness filters. GTDB-Tk will not be executed.")
            }
        }

    GTDBTK_SUMMARY(
        ch_discarded_bins.collect().ifEmpty([]),
        GTDBTK_CLASSIFYWF.out.summary.map { _meta, summary -> summary }.collect(),
        [],
        [],
    )

    ch_summary = channel.empty().mix(GTDBTK_SUMMARY.out.summary)
    ch_versions = ch_versions.mix(GTDBTK_SUMMARY.out.versions)

    emit:
    summary       = ch_summary
    multiqc_files = ch_summary
    versions      = ch_versions
}
