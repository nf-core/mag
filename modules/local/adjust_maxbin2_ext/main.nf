process ADJUST_MAXBIN2_EXT {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*.fa.gz"), emit: renamed_bins

    script:
    """
    if [ -n "${bins}" ]
    then
        for file in ${bins}; do
            [[ \${file} =~ (.*).fasta.gz ]];
            bin="\${BASH_REMATCH[1]}"
            mv \${file} \${bin}.fa.gz
        done
    fi
    """
}
