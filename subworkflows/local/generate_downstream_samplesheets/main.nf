//
// Subworkflow with functionality specific to the nf-core/createtaxdb pipeline
//

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_databases

    main:
    ch_header  = Channel.empty()
    format     = 'csv' // most common format in nf-core
    format_sep = ','

    // TODO --
    // Make your samplesheet channel construct here depending on your downstream
    // pipelines
    if ( params.downstream_pipeline == 'taxprofiler' ) {
        format = 'csv'
        format_sep = ','
        ch_list_for_samplesheet = ch_databases
                                    .map {
                                        meta, db ->
                                            def tool      = meta.tool
                                            def db_name   = meta.id + '-' + meta.tool
                                            def db_params = ""
                                            def db_type   = ""
                                            def db_path   = file(params.outdir).getParent() + '/' + meta.tool + '/' + db.getName()
                                        [ tool: tool, db_name: db_name, db_params: db_params, db_type: db_type, db_path: db_path ]
                                    }
                                    .tap{ ch_header }
    }
    // -- FINISH TODO

    // Constructs the header string and then the strings of each row, and
    // finally concatenates for saving.
    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${params.outdir}/downstream_samplesheet/${params.downstream_pipeline}.${format}",
            newLine: true,
            sort: false
        )

}
