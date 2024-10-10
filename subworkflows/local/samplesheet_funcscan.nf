workflow {
    take:
    ch_assemblies

    main:
    ch_list_for_samplesheet = ch_input
                            .map {
                                meta, filename ->
                                    def sample = meta.id
                                    def fasta  = file(params.outdir).toString() + '/Assembly/' + meta.assembler + '/' + filename.getName()
                                [ sample: sample, fasta: fasta ]
                            }
                            .tap{ ch_header }

    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${params.outdir}/downstream_samplesheet/funcsacn.csv",
            newLine: true,
            sort: false
        )
}
