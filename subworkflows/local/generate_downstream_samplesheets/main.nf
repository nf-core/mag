//
// Subworkflow with functionality specific to the nf-core/mag pipeline
//

workflow SAMPLESHEET_TAXPROFILER {
    take:
    ch_reads

    main:
    format     = 'csv' // most common format in nf-core
    format_sep = ','

    def fastq_rel_path = '/'
    if (params.bbnorm) {
        fastq_rel_path = '/bbmap/bbnorm/'
    } else if (!params.keep_phix) {
        fastq_rel_path = '/QC_shortreads/remove_phix/'
    }
    else if (params.host_fasta) {
        fastq_rel_path = '/QC_shortreads/remove_host/'
    }
    else if (!params.skip_clipping) {
        fastq_rel_path = '/QC_shortreads/fastp/'
    }

    ch_list_for_samplesheet = ch_reads
        .map {
            meta, fastq ->
                def sample              = meta.id
                def run_accession       = meta.id
                def instrument_platform = ""
                def fastq_1             = file(params.outdir).toString() + fastq_rel_path + meta.id + '/' + fastq[0].getName()
                def fastq_2             = file(params.outdir).toString() + fastq_rel_path + meta.id + '/' + fastq[1].getName()
                def fasta               = ""
            [ sample: sample, run_accession: run_accession, instrument_platform: instrument_platform, fastq_1: fastq_1, fastq_2: fastq_2, fasta: fasta ]
        }
        .tap{ ch_colnames }

    channelToSamplesheet(ch_colnames, ch_list_for_samplesheet, 'downstream_samplesheets', 'taxprofiler', format, format_sep)

}

workflow SAMPLESHEET_FUNCSCAN {
    take:
    ch_assemblies

    main:
    format     = 'csv' // most common format in nf-core
    format_sep = ','

    ch_list_for_samplesheet = ch_assemblies
        .map {
            meta, filename ->
                def sample = meta.id
                def fasta  = file(params.outdir).toString() + '/Assembly/' + meta.assembler + '/' + filename.getName()
            [ sample: sample, fasta: fasta ]
        }
        .tap{ ch_colnames }

    channelToSamplesheet(ch_colnames, ch_list_for_samplesheet, 'downstream_samplesheets', 'funcscan', format, format_sep)

}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads
    ch_assemblies

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if ( downstreampipeline_names.contains('taxprofiler') && params.save_clipped_reads ) { // save_clipped_reads must be true
        SAMPLESHEET_TAXPROFILER(ch_reads)
    }

    if ( downstreampipeline_names.contains('funcscan') ) {
        SAMPLESHEET_FUNCSCAN(ch_assemblies)
    }
}

// Constructs the header string and then the strings of each row, and
def channelToSamplesheet(ch_header, ch_list_for_samplesheet, outdir_subdir, pipeline, format, format_sep) {
    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${params.outdir}/${outdir_subdir}/${pipeline}.${format}",
            newLine: true,
            sort: false
        )
}
