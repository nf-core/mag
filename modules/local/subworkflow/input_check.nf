/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include {
    SAMPLESHEET_CHECK;
    get_samplesheet_short_reads_paths;
    get_samplesheet_long_reads_paths } from '../process/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    // TODO what about --input '*_R{1,2}.fastq.gz' ?
    SAMPLESHEET_CHECK ( samplesheet )

    SAMPLESHEET_CHECK.out
        .splitCsv ( header:true, sep:',' )
        .map { row ->
            get_samplesheet_short_reads_paths(row)
        }
        .set { short_reads }

    // TODO could use multiMap {} here, but seems that channels must get same number of items
    // (not always long reads given and not sure if this should be handled differently)
    SAMPLESHEET_CHECK.out
        .splitCsv ( header:true, sep:',' )
        .map { row ->
            if (row.long_reads) get_samplesheet_long_reads_paths(row)
        }
        .set { long_reads }

    // TODO if hybrid: check other cases, output warnings here?

    emit:
    short_reads // channel: [ val(meta), [ reads ] ]
    long_reads
}
