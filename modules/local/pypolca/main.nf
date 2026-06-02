process PYPOLCA {
    tag "${meta.assembler}-${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://quay.io/biocontainers/pypolca:0.4.0--pyhdfd78af_0' : 'quay.io/biocontainers/pypolca:0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}_polished.fasta"), emit: polished
    path "versions.yml", emit: versions

    script:
    """
    pypolca run -a ${assembly} -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} -o pypolca_output --careful

    cp pypolca_output/pypolca_corrected.fasta ${meta.id}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
          pypolca: \$(pypolca --version)
    END_VERSIONS
    """
}

