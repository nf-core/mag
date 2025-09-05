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
    ch_contig_bin_assignments = ch_input_for_postbinning
        .transpose()
        .map { meta, binfile -> [meta + [bin_id: binfile.name], binfile] }
        .splitFasta(record: [header: true])
        .map { meta, record -> [[id: meta.id, contig_id: record.header], [bin_id: meta.bin_id]] }

    // Make useable sample id/contig name meta map for merging
    ch_pydamage_filtered_results = ch_pydamage_results
        .splitCsv(header: true)
        .map { meta, pydamage_stats -> [[id: meta.id, contig_id: pydamage_stats.reference], meta, pydamage_stats] }

    //Merge sample ID to pydamage results so we can group by contig AND node name for 
    ch_reordered_pydamage_stats = ch_contig_bin_assignments
        .join(ch_pydamage_filtered_results)
        .map { _coremeta, _bin_meta, _full_meta, pydamage_stats -> [_bin_meta, pydamage_stats] }
        .groupTuple(by: 0)
        .dump(tag: 'grouped')
        .multiMap { meta, data ->
            bin_id: meta.bin_id
            contents: data
        }

    // TODO: somehow extract bin_id as a string not DataBroardcast thing
    // TODO: check the resulting file HAS reorganised the contigs as expected (per bin)
    ch_pydamage_to_bins = channelToTable(ch_reordered_pydamage_stats.contents, "${params.outdir}/Ancient_DNA/pydamage/bin_summary/${ch_reordered_pydamage_stats.bin_id[0].toString()}", 'csv')

    // TODO: do a dummy test to see if you can pass the collectFile as an object to e.g. a process
    // SUMMARISEPYDAMAGE([[id: 'test'], ch_pydamage_to_bins])
    ch_versions = ch_versions.mix(Channel.empty())

    emit:
    tsv      = Channel.empty() // SUMMARISEPYDAMAGE.out.summary_tsv
    versions = ch_versions
}


// Constructs the header string and then the strings of each row, and
def channelToTable(ch_list_for_samplesheet, path, format) {
    def format_sep = ["csv": ",", "tsv": "\t", "txt": "\t"][format]

    def ch_header = ch_list_for_samplesheet

    def final_file = ch_header
        .first()
        .view()
        .map { it[0].keySet().join(format_sep) }
        .concat(ch_list_for_samplesheet.map { it[0].values().join(format_sep) })
        .view()
        .collectFile(
            name: "${path}.${format}",
            newLine: true,
            sort: false,
        )

    return final_file
}
