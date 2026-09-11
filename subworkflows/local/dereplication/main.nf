//
// Dereplicate bins study-wide with Galah, using each bin's best-available
// completeness/contamination estimate (from whichever of BUSCO/CheckM/CheckM2
// ran) to pick a representative genome per cluster.
//
// Bin filenames are prefixed per-sample/binner upstream, so they're unique
// study-wide and safe to use as the join key back to per-sample metadata.
//

include { GALAH } from '../../../modules/nf-core/galah/main'

workflow DEREPLICATION {

    take:
    ch_bins        // channel: [ val(meta), [ path(bin) ] ], per-sample bins (mandatory)
    ch_genome_info // channel: path(csv), single study-wide "genome,completeness,contamination" table from BIN_QC.out.genome_info (bare path, no meta)

    main:
    // Unbinned-contig pseudo-bins aren't real genomes -- QC tools don't
    // reliably assess them, and Galah panics on any bin missing from its
    // QC report -- so exclude them here rather than crashing downstream.
    ch_bins_flat = ch_bins
        .filter { meta, _bins -> meta.domain != "eukarya" && !meta.refinement.endsWith("unbinned") }
        .transpose()

    ch_bins_for_galah = ch_bins_flat
        .map { _meta, bin -> bin }
        .collect()
        .map { bins -> [[id: 'study'], bins] }

    ch_galah_input = ch_bins_for_galah
        .combine(ch_genome_info)
        .map { meta, bins, qc -> [meta, bins, qc, 'genome-info'] }

    // Galah panics on zero qualifying genomes instead of erroring cleanly
    // (https://github.com/wwood/galah/issues/75), so count them ourselves
    // first and skip Galah gracefully instead.
    ch_qualifying_count = ch_genome_info.map { csv ->
        csv.splitCsv(header: true).count { row ->
            (row.completeness as Double) >= params.dereplicate_min_completeness && (row.contamination as Double) <= params.dereplicate_max_contamination
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

    // Pass every bin through unclustered rather than starving
    // GTDB-Tk/CAT-BAT/Prokka when Galah was skipped.
    ch_fallback_bins = ch_galah_routed.skip
        .combine(ch_bins_flat)
        .map { _count, meta, bin -> [meta, bin] }

    // Recover each representative's original per-sample metadata by joining
    // back on filename.
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
