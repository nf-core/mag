process CONVERT_DEPTHS {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/82/8266909aff235d1f423ca76029d89a7b470f17c634e001aa388b85713d71dd16/data' :
        'community.wave.seqera.io/library/bioawk:1.0--9a4151b4c4cbe47a' }"

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    // need to add empty val because representing reads as we dont want maxbin to calculate for us.
    tuple val(meta), path(fasta), val([]), path("*.abund"), emit: output
    path "versions.yml"                                   , emit: versions

    script:
    """
    # Write one abundance file per read set in a single streaming pass. The depth
    # file is never decompressed to disk: on shared/object-backed work dirs the
    # plaintext copy and one full read per column are prohibitively slow.
    gzip -cd ${depth} | bioawk -t '
        NR == 1 {
            n_abund = int((NF - 3) / 2)
            for (i = 1; i <= n_abund; i++) {
                name = \$(i * 2 + 2)
                sub(/\\.bam\$/, "", name)
                abund[i] = name ".abund"
            }
            next
        }
        {
            for (i = 1; i <= n_abund; i++) {
                print \$1, \$(i * 2 + 2) > (abund[i])
            }
        }
    '

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
