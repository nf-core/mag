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
    tuple val(meta), path("*.tooShort.fa.gz")         , optional:true , emit: tooshort
    tuple val(meta), path("*.lowDepth.fa.gz")         , optional:true , emit: lowdepth
    tuple val(meta), path("*.unbinned.fa.gz")         , optional:true , emit: unbinned
    tuple val(meta), path("*.tsv.gz")                 , optional:true , emit: membership
    tuple val(meta), path("MetaBinner_bins/*.fa.gz")  , optional:true , emit: bins
    tuple val(meta), path("*.log.gz")                 , optional:true , emit: log
    tuple val(meta), path("coverage_profile.tsv")     , optional:true , emit: coverage_profile
    tuple val(meta), path("*kmer_4_f1000.csv")        , optional:true , emit: composition_profile
    path "versions.yml"                               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metabinner_path=\$(dirname \$(which run_metabinner.sh))
    f="${fasta}"
    outname=\${f%.*}
    wd=\$(pwd)

    # create coverage profile in Metabinner format
    zcat \${wd}/${depth} | awk '{if (\$2>${args}) print \$0 }' | cut -f -1,4- > coverage_profile.tsv

    # create composition profile (contigs > ${args} p (default 1000), k = 4)
    python \${metabinner_path}/scripts/gen_kmer.py ${fasta} ${args} 4

    # filter contigs < ${args} bp from fasta (default 1000)
    python \${metabinner_path}/scripts/Filter_tooshort.py ${fasta} ${args}

    bash run_metabinner.sh \\
        -a \${wd}/\${outname}_${args}.fa \\
        -o \${wd}/${prefix} \\
        -d \${wd}/coverage_profile.tsv \\
        -k \${wd}/\${outname}_kmer_4_f${args}.csv \\
        -t ${task.cpus} \\
        -p \${metabinner_path}

    mv ${prefix}/metabinner_res/metabinner_result.tsv ${prefix}.tsv
    python create_metabinner_bins.py ${prefix}.tsv ${fasta} \${wd}/metabinner_bins ${prefix}

    gzip -cn ${prefix} > ${prefix}.tsv.gz

    mv \${outname}_${args}.fa \${outname}.tooShort.fa
    find ./metabinner_bins/ -name "*.fa" -type f | xargs -t -n 1 bgzip -@ ${task.cpus}
    find ./metabinner_bins/ -name "*[lowDepth,tooShort,unbinned].fa.gz" -type f -exec mv {} . \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaBinner: 1.4.4-0
    END_VERSIONS
    """
}
