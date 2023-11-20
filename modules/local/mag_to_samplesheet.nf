process MAG_TO_SAMPLESHEET {
    tag "$meta.id"

    executor 'local'
    memory 100.MB

    input:
    val meta
    val pipeline

    output:
    tuple val(meta), path("*samplesheet.csv"), emit: samplesheet

    exec:
    //
    // Create samplesheet containing metadata
    //

    // Add nf-core pipeline specific entries
    if (pipeline) {
        if (pipeline == 'phageannotator') {
            pipeline_map = [
                sample  : "${meta.id}",
                group   : "${meta.group}",
                fastq_1 : meta.fastq_1,
                fastq_2 : meta.fastq_2,
                fasta   : meta.fasta
            ]
        }
    }

    // Create a samplesheet
    samplesheet  = pipeline_map.keySet().collect{ '"' + it + '"'}.join(",") + '\n'
    samplesheet += pipeline_map.values().collect{ '"' + it + '"'}.join(",")

    // Write samplesheet to file
    def samplesheet_file = task.workDir.resolve("${meta.id}.samplesheet.csv")
    samplesheet_file.text = samplesheet

}
