process MAG_MERGE_SAMPLESHEET {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path ('samplesheets/*')

    output:
    path "*_samplesheet.csv", emit: samplesheet
    path "versions.yml"   , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 1 `ls ./samplesheets/* | head -n 1` > ${prefix}_samplesheet.csv
    for fileid in `ls ./samplesheets/*`; do
        awk 'NR>1' \$fileid >> ${prefix}_samplesheet.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
