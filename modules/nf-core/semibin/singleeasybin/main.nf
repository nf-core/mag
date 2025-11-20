process SEMIBIN_SINGLEEASYBIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/semibin:2.2.0--pyhdfd78af_0':
        'biocontainers/semibin:2.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(bam)

    output:
    tuple val(meta), path("${prefix}/*.csv")              , emit: csv
    tuple val(meta), path("${prefix}/*.h5")               , emit: model, optional: true
    tuple val(meta), path("${prefix}/*.tsv")              , emit: tsv
    tuple val(meta), path("*.log")                        , emit: log
    tuple val(meta), path("bins/${prefix}*.fa.gz")        , emit: output_fasta, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    prefix    = task.ext.prefix ?: "${meta.id}"
    """

    SemiBin2 \\
        $args \\
        single_easy_bin \\
        --input-fasta ${fasta} \\
        --input-bam ${bam} \\
        --tag-output ${prefix} \\
        --output ${prefix} \\
        -t $task.cpus \\
        $args2

    # move final bins to "bins" folder
    mv ./${prefix}/output_bins/ bins/
    # rename log with prefix
    mv ./${prefix}/SemiBinRun.log ${prefix}.SemiBinRun.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin: \$( SemiBin2 --version )
    END_VERSIONS
    """
    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/{contig_bins,recluster_bins_info}.tsv
    touch ${prefix}/{data,data_split}.csv
    mkdir bins
    touch bins/${prefix}_SemiBin_{0,1,2,3}.fa
    gzip  bins/${prefix}_SemiBin*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin: \$( SemiBin2 --version )
    END_VERSIONS
    """
}
