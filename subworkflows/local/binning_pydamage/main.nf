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
    // 5. Final emitted TSV is then used in mag.nf bind to bin_summary table
    //     /*

    ch_collected_pydamage_results = ch_contig_pydamage_results
        .collect { _meta, pydamage_report -> pydamage_report }
        .map { pydamage_reports -> pydamage_reports.sort() }

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
        )
        .dump(tag: 'collected_file')

    SUMMARISE_PYDAMAGEBINS(ch_collected_pydamage_results, ch_bin_contig_names)

    emit:
    tsv      = channel.empty()
    versions = ch_versions
}
