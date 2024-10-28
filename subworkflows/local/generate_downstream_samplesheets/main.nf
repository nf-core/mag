//
// Subworkflow with functionality specific to the nf-core/mag pipeline
//

workflow SAMPLESHEET_TAXPROFILER {
    take:
    ch_reads

    main:
    format = 'csv'

    def fastq_rel_path = '/'
    if (params.bbnorm) {
        fastq_rel_path = "/bbmap/bbnorm/"
    }
    else if (!params.keep_phix) {
        fastq_rel_path = "/QC_shortreads/remove_phix/"
    }
    else if (params.host_fasta != false) {
        fastq_rel_path = "/QC_shortreads/remove_host/"
    }
    else if (!params.skip_clipping && params.clip_tool == 'fastp') {
        fastq_rel_path = "/QC_shortreads/fastp/"
    }
    else if (!params.skip_clipping && params.clip_tool == 'adapterremoval') {
        fastq_rel_path = "/QC_shortreads/adapterremoval/"
    }

    ch_list_for_samplesheet = ch_reads
        .map { meta, fastq ->
            def sample = meta.id
            def run_accession = meta.id
            def instrument_platform = ""
            def fastq_1 = meta.single_end ? file(params.outdir).toString() + fastq_rel_path + meta.id + '/' + fastq.getName() : file(params.outdir).toString() + fastq_rel_path + meta.id + '/' + fastq[0].getName()
            def fastq_2 = meta.single_end ? "" : file(params.outdir).toString() + fastq_rel_path + meta.id + '/' + fastq[1].getName()
            def fasta = ""
            [sample: sample, run_accession: run_accession, instrument_platform: instrument_platform, fastq_1: fastq_1, fastq_2: fastq_2, fasta: fasta]
        }
        .tap { ch_colnames }

    channelToSamplesheet(ch_list_for_samplesheet, "${params.outdir}/downstream_samplesheets/taxprofiler", format)
}

workflow SAMPLESHEET_FUNCSCAN {
    take:
    ch_assemblies

    main:
    format = 'csv'

    ch_list_for_samplesheet = ch_assemblies
        .map { meta, filename ->
            // funcscan requires
            def sample = filename.extension ==~ 'gz' ? filename.baseName.take(filename.baseName.lastIndexOf('.')) : filename.baseName
            def fasta = file(params.outdir).toString() + '/Assembly/' + meta.assembler + '/' + filename.getName()
            [sample: sample, fasta: fasta]
        }
        .tap { ch_colnames }

    channelToSamplesheet(ch_list_for_samplesheet, "${params.outdir}/downstream_samplesheets/funcscan", format)
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads
    ch_assemblies

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if (downstreampipeline_names.contains('taxprofiler')) {
        SAMPLESHEET_TAXPROFILER(ch_reads)
    }

    if (downstreampipeline_names.contains('funcscan')) {
        SAMPLESHEET_FUNCSCAN(ch_assemblies)
    }
}

// Constructs the header string and then the strings of each row, and
def channelToSamplesheet(ch_list_for_samplesheet, path, format) {
    def format_sep = [csv: ",", tsv: "\t", txt: "\t"][format]

    def ch_header = ch_list_for_samplesheet

    ch_header
        .first()
        .map { it.keySet().join(format_sep) }
        .concat(ch_list_for_samplesheet.map { it.values().join(format_sep) })
        .collectFile(
            name: "${path}.${format}",
            newLine: true,
            sort: false
        )
}
