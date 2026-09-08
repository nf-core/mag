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

        GALAH(ch_galah_input)

        ch_representative_files = GALAH.out.dereplicated_bins
        ch_cluster_tsv = GALAH.out.tsv
    }

    // Recover each representative's original per-sample metadata by
    // joining back on filename against the pre-collection flat channel.
    ch_bins_keyed = ch_bins_flat.map { meta, bin -> [bin.name, meta, bin] }

    ch_representative_keys = ch_representative_files
        .transpose()
        .map { _meta, bin -> [bin.name, bin] }

    ch_dereplicated_bins = ch_representative_keys
        .join(ch_bins_keyed)
        .map { _name, _representative_bin, meta, bin -> [meta, bin] }

    // Galah emits its version only via the pipeline-wide versions topic
    // channel, so there's no per-subworkflow versions channel to emit here.
    emit:
    dereplicated_bins = ch_dereplicated_bins // channel: [ val(meta), path(bin) ], one representative genome per cluster, original metadata preserved
    cluster_tsv       = ch_cluster_tsv // channel: [ val(meta), path(tsv) ], representative <TAB> member cluster definition
}
