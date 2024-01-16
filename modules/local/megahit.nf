process MEGAHIT {
    tag "$meta.id"

    conda "bioconda::megahit=1.2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megahit:1.2.9--h2e03b76_1' :
        'biocontainers/megahit:1.2.9--h2e03b76_1' }"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("MEGAHIT/MEGAHIT-${meta.id}.contigs.fa"), emit: assembly
    path "MEGAHIT/*.log"                                  , emit: log
    path "MEGAHIT/MEGAHIT-${meta.id}.contigs.fa.gz"               , emit: assembly_gz
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def input = meta.single_end ? "-r \"" + reads1.join(",") + "\"" : "-1 \"" + reads1.join(",") + "\" -2 \"" + reads2.join(",") + "\""
    mem = task.memory.toBytes()
    if ( !params.megahit_fix_cpu_1 || task.cpus == 1 )
        """
        ## Check if we're in the same work directory as a previous failed MEGAHIT run
        if [[ -d MEGAHIT ]]; then
            rm -r MEGAHIT/
        fi

        megahit $args -t "${task.cpus}" -m $mem $input -o MEGAHIT --out-prefix "MEGAHIT-${meta.id}"

        gzip -c "MEGAHIT/MEGAHIT-${meta.id}.contigs.fa" > "MEGAHIT/MEGAHIT-${meta.id}.contigs.fa.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    else
        error "ERROR: '--megahit_fix_cpu_1' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}
