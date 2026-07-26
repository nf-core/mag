process TIARA_CLASSIFY {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1021c2bc41756fa99bc402f461dad0d1c35358c1:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' :
        'biocontainers/mulled-v2-1021c2bc41756fa99bc402f461dad0d1c35358c1:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' }"

    input:
    tuple val(meta), path(classification), path(contig2bin), path(bins)

    output:
    tuple val(meta), path("eukarya/*.fa.gz"),         emit: eukarya_bins, optional: true
    tuple val(meta), path("prokarya/*.fa.gz"),        emit: prokarya_bins, optional: true
    tuple val(meta), path("bacteria/*.fa.gz"),        emit: bacteria_bins, optional: true
    tuple val(meta), path("archaea/*.fa.gz"),         emit: archaea_bins, optional: true
    tuple val(meta), path("organelle/*.fa.gz"),       emit: organelle_bins, optional: true
    tuple val(meta), path("unknown/*.fa.gz"),         emit: unknown_bins, optional: true
    tuple val(meta), path("*.binclassification.tsv"), emit: bin_classifications
    path 'versions.yml',                              emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    domain_classification.R \
        --classification_file ${classification} \
        --contig_to_bin ${contig2bin} \
        ${args} \
        --output_prefix ${prefix}

    mkdir eukarya
    mkdir prokarya
    mkdir bacteria
    mkdir archaea
    mkdir organelle
    mkdir unknown

    while IFS=\$"\t" read bin domain; do
        # dereference with -L: staged inputs are symlinks, and moving one would
        # emit a dangling link instead of the bin contents on remote filesystems.
        # match the extension exactly, otherwise e.g. bin .1 also picks up .10
        find -L . -maxdepth 1 -name "\${bin}.fa.gz" -exec cp -L {} \${domain}/ \\;
    done < bin2classification.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | grep -Eo '[0-9.]+ ')
        r-tidyverse: \$(cat tidyverse_version.txt)
    END_VERSIONS
    """
}
