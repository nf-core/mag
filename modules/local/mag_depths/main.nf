process MAG_DEPTHS {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"

    // Using container from metabat2 process, since this will be anyway already downloaded and contains biopython and pandas
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' :
        'biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' }"

    input:
    tuple val(meta), path(bins), path(contig_depths)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.binner}-${meta.id}-binDepths.tsv"), emit: depths
    path "versions.yml"                                                               , emit: versions

    script:
    """
    get_mag_depths.py --bins ${bins} \\
                    --depths ${contig_depths} \\
                    --assembler ${meta.assembler} \\
                    --id ${meta.id} \\
                    --binner ${meta.binner}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
