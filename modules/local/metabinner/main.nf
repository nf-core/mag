process METABINNER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabinner:1.4.4--hdfd78af_0' :
        'quay.io/biocontainers/metabinner:1.4.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    tuple val(meta), path("*.tooShort.fa.gz")         , emit: tooshort
    tuple val(meta), path("*.lowDepth.fa.gz")         , optional:true , emit: lowdepth
    tuple val(meta), path("*.unbinned.fa.gz")         , emit: unbinned
    tuple val(meta), path("*.tsv.gz")                 , emit: membership
    tuple val(meta), path("metabinner_bins/*.fa.gz")  , emit: bins
    tuple val(meta), path("*.log.gz")                 , emit: log
    tuple val(meta), path("coverage_profile.tsv")     , emit: coverage_profile
    tuple val(meta), path("*kmer_4_f1000.csv")        , emit: composition_profile
    path "versions.yml"                               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_contig_size = task.ext.min_contig_size ?: "1000"
    """
    metabinner_path=\$(dirname \$(which run_metabinner.sh))
    f="${fasta}"

    # create coverage profile in Metabinner format
    zcat ${depth} | awk '{if (\$2>${min_contig_size}) print \$0 }' | cut -f -1,4- > coverage_profile.tsv

    # create composition profile (contigs > ${min_contig_size} p (default 1000), k = 4)
    python \${metabinner_path}/scripts/gen_kmer.py ${fasta} ${min_contig_size} 4

    # filter contigs < ${min_contig_size} bp from fasta (default 1000)
    python \${metabinner_path}/scripts/Filter_tooshort.py ${fasta} ${min_contig_size}

    # requires absolute paths
    wd=\$(pwd)
    outname=\${f%.*}
    bash run_metabinner.sh \\
        -a \${wd}/\${outname}_${min_contig_size}.fa \\
        -o \${wd}/${prefix} \\
        -d \${wd}/coverage_profile.tsv \\
        -k \${wd}/\${outname}_kmer_4_f${min_contig_size}.csv \\
        -t ${task.cpus} \\
        -p \${metabinner_path} \\
        ${args}

    # collect & zip membership & log files
    gzip -cn ${prefix}/metabinner_res/metabinner_result.tsv > ${prefix}.tsv.gz
    gzip -cn ${prefix}/metabinner_res/result.log > ${prefix}.log.gz

    # collect & zip bins & un-binned fractions
    create_metabinner_bins.py \\
        ${prefix}/metabinner_res/metabinner_result.tsv \\
        ${fasta} \\
        ./metabinner_bins \\
        ${prefix} \\
        ${min_contig_size}
    find ./metabinner_bins/ -name "*.fa" -type f | xargs -t -n 1 bgzip -@ ${task.cpus}

    # collect & zip non-binned fractions
    find ./metabinner_bins/ -name "*[lowDepth,tooShort,unbinned].fa.gz" -type f -exec mv {} . \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaBinner: 1.4.4-0
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
