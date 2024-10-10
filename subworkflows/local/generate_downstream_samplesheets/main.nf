//
// Subworkflow with functionality specific to the nf-core/createtaxdb pipeline
//

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads
    ch_assemblies
    downstreampipeline_name

    main:
    format     = 'csv' // most common format in nf-core
    format_sep = ','
    // Make your samplesheet channel construct here depending on your downstream pipelines
    if ( downstreampipeline_name == 'taxprofiler' && params.save_clipped_reads ) { // save_clipped_reads must be true
        SMAPLESHEET_TAXPROFILR(ch_reads)
    }

    if ( downstreampipeline_name == 'funcscan' ) {
        SAMPLESHEET_FUNCSACN(ch_assemblies)
    }

}
