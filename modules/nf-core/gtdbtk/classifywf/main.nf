process GTDBTK_CLASSIFYWF {
    tag "${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::gtdbtk=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.1.1--pyhdfd78af_1' :
        'biocontainers/gtdbtk:2.1.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*")

    output:
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.*.summary.tsv"        , emit: summary
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.*.classify.tree.gz"   , emit: tree
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.*.markers_summary.tsv", emit: markers
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.*.msa.fasta.gz"       , emit: msa
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.*.user_msa.fasta.gz"  , emit: user_msa
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.*.filtered.tsv"       , emit: filtered
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.log"                  , emit: log
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.warnings.log"         , emit: warnings
    path "gtdbtk.${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}.failed_genomes.tsv"   , emit: failed
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pplacer_scratch = params.gtdbtk_pplacer_scratch ? "--scratch_dir pplacer_tmp" : ""
    def prefix = task.ext.prefix ?: "${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}"

    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    if [ ${pplacer_scratch} != "" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf \\
        $args \\
        --genome_dir bins \\
        --prefix "gtdbtk.${prefix}" \\
        --out_dir "\${PWD}" \\
        --cpus $task.cpus \\
        --pplacer_cpus $params.gtdbtk_pplacer_cpus \\
        $pplacer_scratch \\
        --min_perc_aa $params.gtdbtk_min_perc_aa \\
        --min_af $params.gtdbtk_min_af

    mv classify/gtdbtk.${prefix}.*.classify.tree \\
        classify/gtdbtk.${prefix}.*.summary.tsv \\
        .
    
    mv identify/gtdbtk.${prefix}.*.markers_summary.tsv \\
       identify/gtdbtk.${prefix}.failed_genomes.tsv \\
       .

    mv align/gtdbtk.${prefix}.*.msa.fasta.gz \\
       align/gtdbtk.${prefix}.*.user_msa.fasta.gz \\
       align/gtdbtk.${prefix}.*.filtered.tsv \\
       . 
       
    gzip gtdbtk.${prefix}.*.classify.tree

    

    mv gtdbtk.log "gtdbtk.${prefix}.log"
    mv gtdbtk.warnings.log "gtdbtk.${prefix}.warnings.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    def VERSION = '2.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch gtdbtk.${prefix}.stub.summary.tsv
    touch gtdbtk.${prefix}.stub.classify.tree.gz
    touch gtdbtk.${prefix}.stub.markers_summary.tsv
    touch gtdbtk.${prefix}.stub.msa.fasta.gz
    touch gtdbtk.${prefix}.stub.user_msa.fasta.gz
    touch gtdbtk.${prefix}.stub.filtered.tsv
    touch gtdbtk.${prefix}.log
    touch gtdbtk.${prefix}.warnings.log
    touch gtdbtk.${prefix}.failed_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo "$VERSION")
    END_VERSIONS
    """
}
