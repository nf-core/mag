process CONVERT_DEPTHS {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--hed695b0_5' :
        'biocontainers/bioawk:1.0--hed695b0_5' }"

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    // need to add empty val because representing reads as we dont want maxbin to calculate for us.
    tuple val(meta), path(fasta), val([]), path("*.abund"), emit: output
    path "versions.yml"                                   , emit: versions

    script:
    """
    gunzip -f ${depth}

    # Determine the number of abundance columns
    n_abund=\$(awk 'NR==1 {print int((NF-3)/2)}' ${depth.toString() - '.gz'})

    # Get column names
    read -r header<${depth.toString() - '.gz'}
    header=(\$header)

    # Generate abundance files for each read set
    for i in \$(seq 1 \$n_abund); do
        col=\$((i*2+2))
        name=\$( echo \${header[\$col-1]} | sed s/\\.bam\$// )
        bioawk -t '{if (NR > 1) {print \$1, \$'"\$col"'}}' ${depth.toString() - '.gz'} > \${name}.abund
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
