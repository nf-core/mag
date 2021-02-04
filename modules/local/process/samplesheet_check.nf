// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet
    
    output:
    path '*.csv'


    script:  // This script is bundled with the pipeline, in nf-core/mag/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ meta, [ sr1, sr2 ] ]
def get_samplesheet_short_reads_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.group        = row.group
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.short_reads_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.short_reads_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.short_reads_1) ] ]
    } else {
        if (!file(row.short_reads_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.short_reads_2}"
        }
        array = [ meta, [ file(row.short_reads_1), file(row.short_reads_2) ] ]
    }
    return array
}

// Function to get list of [ meta, [ lr ] ]
def get_samplesheet_long_reads_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.group        = row.group

    def array = []
    if (row.long_reads && !file(row.long_reads).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Long Read FastQ file does not exist!\n${row.long_reads}"
    } else if (row.long_reads) {
        array = [ meta, [ file(row.long_reads) ] ]
    }
    return array
}
