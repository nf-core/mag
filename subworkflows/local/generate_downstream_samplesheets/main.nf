//
// Subworkflow with functionality specific to the nf-core/createtaxdb pipeline
//

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads

    main:
    format     = 'csv' // most common format in nf-core
    format_sep = ','
    // Make your samplesheet channel construct here depending on your downstream pipelines
    if ( params.generate_pipeline_samplesheets == 'taxprofiler' && params.save_clipped_reads ) { // save_clipped_reads must be true
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
            .tap{ ch_header }
    }

    // Constructs the header string and then the strings of each row, and
    // finally concatenates for saving.
    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${params.outdir}/downstream_samplesheet/${params.generate_pipeline_samplesheets}.${format}",
            newLine: true,
            sort: false
        )

}
