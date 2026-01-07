include { SUMMARISE_PYDAMAGEBINS } from '../../../modules/local/summarise_pydamagebins/main'

workflow BINNING_PYDAMAGE {
    take:
    ch_contig_pydamage_results
    ch_input_for_postbinning

    main:
    ch_versions = channel.empty()

    // 1. Collect all pydamage results into single tuple
    // 2. Generate list of contigs per bin
    // 3. Local module that reorders/reassigns contigs (save), and also summarises pydamage results per bin via median
    // 4. Emit final tsv and versions
    // DOWNSTREAM: Final emitted TSV is then used in `mag.nf` to bind to final `bin_summary.tsv` table

    ch_collected_pydamage_results = ch_contig_pydamage_results
        .map { _meta, pydamage_report -> pydamage_report }
        .toSortedList { pydamage_report ->
            // Sort based on filename only (not full path) as work directory path will be different each run
            file(pydamage_report).getBaseName()
        }

    ch_bin_contig_names = ch_input_for_postbinning
        .transpose()
        .map { meta, binfile -> [meta + [bin_id: binfile.name], binfile] }
        .splitFasta(record: [header: true], elem: 1)
        .map { meta, contig_header ->
            "${meta['bin_id']}\t${meta['assembler']}-${meta['id']}\t${meta['binner']}\t${contig_header['header']}"
        }
        .collectFile(
            name: 'contig_to_bin_map.tsv',
            storeDir: params.outdir + '/GenomeBinning/QC',
            newLine: true,
            sort: true,
        )

    SUMMARISE_PYDAMAGEBINS(ch_collected_pydamage_results, ch_bin_contig_names)
    ch_versions = ch_versions.mix(SUMMARISE_PYDAMAGEBINS.out.versions)

    emit:
    tsv      = SUMMARISE_PYDAMAGEBINS.out.pydamage_bin_summary
    versions = ch_versions
}
