//
// Dereplicate bins study-wide, using CheckM2 quality estimates to pick a
// representative genome per cluster.
//
// Bin filenames are prefixed per-sample/binner upstream (e.g.
// "${meta.assembler}-MetaBAT2-${meta.id}"), so they stay unique study-wide
// and can safely be used as the join key to recover each representative's
// original per-sample metadata after the study-wide clustering call
// collapses everything into a single, unkeyed collection.
//

include { GAWK  } from '../../../modules/nf-core/gawk/main'
include { GALAH } from '../../../modules/nf-core/galah/main'

workflow DEREPLICATION {

    take:
    ch_bins            // channel: [ val(meta), [ path(bin) ] ], per-sample bins (mandatory)
    ch_checkm2_summary // channel: path(tsv), single study-wide CheckM2 summary from BIN_QC.out.checkm2_summary (bare path, no meta)

    main:
    // Leftover unbinned-contig pseudo-bins (meta.refinement ends in
    // "unbinned": 'unrefined_unbinned' or 'dastool_refined_unbinned') flow
    // into ch_input_for_postbinning alongside real curated bins, but aren't
    // real genomes -- CheckM2 doesn't reliably produce QC for them, and
    // Galah panics outright (rather than warning) on any bin missing from
    // its QC report, so they're excluded here rather than just being an
    // occasional dereplication no-op.
    ch_bins_flat = ch_bins
        .filter { meta, _bins -> meta.domain != "eukarya" && !meta.refinement.endsWith("unbinned") }
        .transpose()

    // Declared here rather than inside the `if` below for the same reason
    // documented further down for GALAH.out.* -- a variable first assigned
    // inside an `if` isn't reliably visible outside it under Nextflow's
    // strict-syntax parser.
    ch_fallback_bins = channel.empty()

    if (params.dereplication_tool == "galah") {
        ch_bins_for_galah = ch_bins_flat
            .map { _meta, bin -> bin }
            .collect()
            .map { bins -> [[id: 'study'], bins] }

        // CheckM2's real Name column is the bare bin stem with no extension
        // at all, but Galah looks entries up by the bin's original extension
        // (with any .gz compression suffix stripped) -- the same mismatch
        // nf-core/modules' own galah module test works around with a GAWK
        // rewrite step, and for the same reason.
        GAWK(
            ch_checkm2_summary.map { tsv -> [[id: 'checkm2'], tsv] },
            [],
            false,
        )
        ch_checkm2_for_galah = GAWK.out.output.map { _meta, tsv -> tsv }

        ch_galah_input = ch_bins_for_galah
            .combine(ch_checkm2_for_galah)
            .map { meta, bins, qc -> [meta, bins, qc, 'checkm2'] }

        // Galah panics (index out of bounds) instead of erroring cleanly or
        // just emitting an empty result when zero genomes pass its quality
        // thresholds -- reported upstream: https://github.com/wwood/galah/issues/75.
        // Count qualifying genomes ourselves first so a study where every
        // bin happens to fail the threshold doesn't crash the whole
        // pipeline; skip Galah gracefully instead, with a clear warning.
        ch_qualifying_count = ch_checkm2_for_galah.map { tsv ->
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

        // When Galah is skipped above, GALAH.out is never populated, so
        // treat dereplication as a no-op for this run instead of silently
        // starving GTDB-Tk/CAT-BAT/Prokka of any input: every bin becomes
        // its own "representative". ch_galah_routed.skip carries exactly
        // one item (the study-wide qualifying count) only when the skip
        // branch fired, so this combine() broadcasts it across every bin
        // and produces nothing at all when Galah actually ran instead.
        ch_fallback_bins = ch_galah_routed.skip
            .combine(ch_bins_flat)
            .map { _count, meta, bin -> [meta, bin] }
    }

    // Recover each representative's original per-sample metadata by
    // joining back on filename against the pre-collection flat channel.
    // Referencing GALAH.out.* directly here (rather than assigning it to an
    // intermediate variable inside the `if` block above) is deliberate: a
    // variable first assigned inside an `if` isn't reliably visible outside
    // it under Nextflow's strict-syntax parser, even though it always would
    // be at runtime while "galah" is the only dereplication_tool value --
    // same pattern DOMAIN_CLASSIFICATION already uses for TIARA.out.* .
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
