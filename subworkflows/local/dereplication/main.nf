//
// Dereplicate bins study-wide, using CheckM or CheckM2 quality estimates to
// pick a representative genome per cluster.
//
// Bin filenames are prefixed per-sample/binner upstream, so they're unique
// study-wide and safe to use as the join key back to per-sample metadata.
//

include { GAWK  } from '../../../modules/nf-core/gawk/main'
include { GALAH } from '../../../modules/nf-core/galah/main'

workflow DEREPLICATION {

    take:
    ch_bins        // channel: [ val(meta), [ path(bin) ] ], per-sample bins (mandatory)
    ch_qc_summary  // channel: path(tsv), single study-wide CheckM or CheckM2 summary from BIN_QC.out.checkm_summary/checkm2_summary (bare path, no meta)
    val_qc_format  // val(string): 'checkm' or 'checkm2', which format ch_qc_summary is in

    main:
    // Unbinned-contig pseudo-bins aren't real genomes -- CheckM/CheckM2
    // don't reliably QC them, and Galah panics on any bin missing from its
    // QC report -- so exclude them here rather than crashing downstream.
    ch_bins_flat = ch_bins
        .filter { meta, _bins -> meta.domain != "eukarya" && !meta.refinement.endsWith("unbinned") }
        .transpose()

    // Declared outside the `if` for the same Nextflow strict-syntax parser
    // scoping reason as GALAH.out.* below.
    ch_fallback_bins = channel.empty()

    if (params.dereplication_tool == "galah") {
        ch_bins_for_galah = ch_bins_flat
            .map { _meta, bin -> bin }
            .collect()
            .map { bins -> [[id: 'study'], bins] }

        // Both CheckM's (--tab_table) and CheckM2's report have the bin ID
        // (extension-less) in column 1, and Galah looks entries up by the
        // bin's original extension (.gz stripped) -- same fix nf-core/
        // modules' own galah module test uses for the same mismatch.
        GAWK(
            ch_qc_summary.map { tsv -> [[id: val_qc_format], tsv] },
            [],
            false,
        )
        ch_qc_summary_for_galah = GAWK.out.output.map { _meta, tsv -> tsv }

        ch_galah_input = ch_bins_for_galah
            .combine(ch_qc_summary_for_galah)
            .map { meta, bins, qc -> [meta, bins, qc, val_qc_format] }

        // Galah panics on zero qualifying genomes instead of erroring
        // cleanly (https://github.com/wwood/galah/issues/75), so count them
        // ourselves first and skip Galah gracefully instead. CheckM's
        // --tab_table and CheckM2's report both name these columns
        // identically ("Completeness"/"Contamination"), so this reads
        // correctly regardless of val_qc_format.
        ch_qualifying_count = ch_qc_summary_for_galah.map { tsv ->
            tsv.splitCsv(header: true, sep: '\t').count { row ->
                (row.Completeness as Double) >= params.dereplicate_min_completeness && (row.Contamination as Double) <= params.dereplicate_max_contamination
            }
        }

        ch_galah_routed = ch_galah_input
            .combine(ch_qualifying_count)
            .branch { meta, bins, qc, format, count ->
                cluster: count > 0
                    return [meta, bins, qc, format]
                skip: true
                    return count
            }

        ch_galah_routed.skip.subscribe {
            log.warn("[nf-core/mag] Dereplication: no bins passed --dereplicate_min_completeness ${params.dereplicate_min_completeness} / --dereplicate_max_contamination ${params.dereplicate_max_contamination}; skipping Galah for this run (works around a Galah crash on zero qualifying genomes, see https://github.com/wwood/galah/issues/75). Every bin is passed through downstream unclustered, as if --dereplicate had not been set, rather than dropped.")
        }

        GALAH(ch_galah_routed.cluster)

        // If Galah was skipped, GALAH.out stays empty -- pass every bin
        // through unclustered instead of starving GTDB-Tk/CAT-BAT/Prokka.
        // ch_galah_routed.skip only carries an item when skip fired, so
        // this combine() broadcasts across every bin only in that case.
        ch_fallback_bins = ch_galah_routed.skip
            .combine(ch_bins_flat)
            .map { _count, meta, bin -> [meta, bin] }
    }

    // Recover each representative's original per-sample metadata by joining
    // back on filename. GALAH.out.* is referenced directly here rather than
    // via an intermediate variable assigned inside the `if` above: Nextflow's
    // strict-syntax parser doesn't reliably see such a variable outside the
    // `if`, even though it always would be at runtime -- same pattern
    // DOMAIN_CLASSIFICATION uses for TIARA.out.* .
    ch_bins_keyed = ch_bins_flat.map { meta, bin -> [bin.name, meta, bin] }

    ch_representative_keys = GALAH.out.dereplicated_bins
        .transpose()
        .map { _meta, bin -> [bin.name, bin] }

    ch_dereplicated_bins = ch_representative_keys
        .join(ch_bins_keyed)
        .map { _name, _representative_bin, meta, bin -> [meta, bin] }
        .mix(ch_fallback_bins)

    // Galah emits its version only via the pipeline-wide versions topic
    // channel, so there's no per-subworkflow versions channel to emit here.
    emit:
    dereplicated_bins = ch_dereplicated_bins // channel: [ val(meta), path(bin) ], one representative genome per cluster, original metadata preserved
    cluster_tsv       = GALAH.out.tsv // channel: [ val(meta), path(tsv) ], representative <TAB> member cluster definition
}
