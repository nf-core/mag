/*
 * SHORTREAD_PREPROCESSING: Preprocessing and QC for short reads
 */

include { SUMMARISEPYDAMAGE } from '../../../modules/local/summarisepydamage/main'

workflow PYDAMAGE_BINS {
    take:
    ch_pydamage_results
    ch_input_for_postbinning

    main:
    ch_versions = Channel.empty()

    // Get sample ID for each contig (record header)
    ch_contig_bin_assignments = ch_input_for_postbinning.transpose()

    // Get meta only version of input bins for later re-binding
    ch_input_for_bins_metaonly = ch_contig_bin_assignments.map { meta, binfile -> [[bin_id: binfile.name], meta] }

    // Extract contig ID per bin
    ch_contig_bin_assignments_meta = ch_contig_bin_assignments
        .map { meta, binfile -> [meta + [bin_id: binfile.name], binfile] }
        .splitFasta(record: [header: true])
        .map { meta, record -> [[id: meta.id, contig_id: record.header], [bin_id: meta.bin_id]] }

    // Make useable sample id/contig name meta map for merging
    ch_pydamage_filtered_results = ch_pydamage_results
        .splitCsv(header: true)
        .map { meta, pydamage_stats -> [[id: meta.id, contig_id: pydamage_stats.reference], meta, pydamage_stats] }

    //Merge sample ID to pydamage results so we can group by contig AND node name
    ch_reordered_pydamage_stats = ch_contig_bin_assignments_meta
        .join(ch_pydamage_filtered_results)
        .map { _coremeta, _bin_meta, _full_meta, pydamage_stats -> [_bin_meta, pydamage_stats] }
        .groupTuple(by: 0)

    // Convert contents of the reordered contigs to a CSV file, and re-attach meta based on bin_id (i.e., from bin file name)
    ch_pydamage_to_bins = ch_reordered_pydamage_stats
        .map { bin_id, data ->
            // Process your data to create CSV content
            def header = data[0].keySet().join(',')
            def rows = data.collect { it.values().join(',') }.join('\n')
            def content = header + '\n' + rows
            [bin_id, content]
        }
        .collectFile(
            [storeDir: "${params.outdir}/GenomeBinning/QC/pydamage/analyze_bins/"]
        ) { meta, content ->
            ["${meta.bin_id}_pydamagebins.csv", content]
        }
        .map { file ->
            def meta = [:]
            meta.bin_id = file.getName() - '_pydamagebins.csv'
            [meta, file]
        }
        .join(ch_input_for_bins_metaonly)
        .map { bin_id, file, meta ->
            def meta_new = meta + bin_id
            [meta_new, file]
        }

    // Generate the per bin-summary
    SUMMARISEPYDAMAGE(ch_pydamage_to_bins)
    ch_versions = ch_versions.mix(SUMMARISEPYDAMAGE.out.versions)

    // Stick per bin-summary into a single summary CSV
    ch_aggregate_summaries = SUMMARISEPYDAMAGE.out.summary_tsv
        .map { _meta, file -> [file] }
        .splitCsv(header: true, sep: '\t')
        .map { content ->
            content[0].id = content[0].id - '_pydamagebins'
            [[id: 'all'], content.flatten()]
        }
        .groupTuple()
        .map { _meta, contents ->
            // Get keys of first row to act as header
            def header = contents[0][0].keySet().join('\t')
            // Get values of remaining rows for cells
            def rows = contents.collect { it[0].values().join('\t') }.join('\n')
            header + '\n' + rows
        }
        .collectFile(
            name: "pydamage_bin_summary.tsv",
            newLine: true,
            sort: false,
            storeDir: "${params.outdir}/GenomeBinning/QC/",
        )

    emit:
    tsv      = ch_aggregate_summaries
    versions = ch_versions
}
