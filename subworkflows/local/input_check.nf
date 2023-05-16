//
// Check input samplesheet and get read channels
//

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow INPUT_CHECK {
    main:
    if(hasExtension(params.input, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        ch_input_rows = Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() >= 5) {
                        def id = row.sample
                        def run = row.run
                        def group = row.group
                        def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : false
                        def sr2 = row.short_reads_2 ? file(row.short_reads_2, checkIfExists: true) : false
                        def lr = row.long_reads ? file(row.long_reads, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (run != null && run == "") exit 1, "ERROR: Please check input samplesheet -> Column 'run' contains an empty field. Either remove column 'run' or fill each field with a value."
                        if (!sr1) exit 1, "Invalid input samplesheet: short_reads_1 can not be empty."
                        if (!sr2 && lr) exit 1, "Invalid input samplesheet: invalid combination of single-end short reads and long reads provided! SPAdes does not support single-end data and thus hybrid assembly cannot be performed."
                        if (!sr2 && !params.single_end) exit 1, "Invalid input samplesheet: single-end short reads provided, but command line parameter `--single_end` is false. Note that either only single-end or only paired-end reads must provided."
                        if (sr2 && params.single_end) exit 1, "Invalid input samplesheet: paired-end short reads provided, but command line parameter `--single_end` is true. Note that either only single-end or only paired-end reads must provided."
                        return [ id, run, group, sr1, sr2, lr ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects at least 5."
                    }
                }
        // separate short and long reads
        ch_raw_short_reads = ch_input_rows
            .map { id, run, group, sr1, sr2, lr ->
                        def meta = [:]
                        meta.id           = id
                        meta.run          = run == null ? "0" : run
                        meta.group        = group
                        meta.single_end   = params.single_end
                        if (params.single_end)
                            return [ meta, [ sr1] ]
                        else
                            return [ meta, [ sr1, sr2 ] ]
                }
        ch_raw_long_reads = ch_input_rows
            .map { id, run, group, sr1, sr2, lr ->
                        if (lr) {
                            def meta = [:]
                            meta.id           = id
                            meta.run          = run == null ? "0" : run
                            meta.group        = group
                            return [ meta, lr ]
                        }
                }
    } else {
        ch_raw_short_reads = Channel
            .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def meta = [:]
                        meta.id           = row[0]
                        meta.run          = 0
                        meta.group        = 0
                        meta.single_end   = params.single_end
                        return [ meta, row[1] ]
                }
        ch_input_rows = Channel.empty()
        ch_raw_long_reads = Channel.empty()
    }

    // Ensure run IDs are unique within samples (also prevents duplicated sample names)
    // Note: do not need to check for PE/SE mismatch, as checks above do not allow mixing
    ch_input_rows
        .groupTuple(by: 0)
        .map { id, run, group, sr1, sr2, lr -> if( run.size() != run.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample or run IDs (within a sample)! Check samplesheet for sample id: ${id}" } }

    emit:
    raw_short_reads = ch_raw_short_reads
    raw_long_reads  = ch_raw_long_reads
}
