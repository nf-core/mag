// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

/*
 * Bowtie2 for read removal
 */
process BOWTIE2_REMOVAL_ALIGN {
    tag "${meta.id}-${options.suffix}"
    publishDir "${params.outdir}/",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1"
    } else {
        container "quay.io/biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1"
    }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.unmapped*.fastq.gz") , emit: reads
    path  "*.mapped*.read_ids.txt", optional:true , emit: read_ids
    tuple val(meta), path("*.bowtie2.log")        , emit: log
    path  '*.version.txt'                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    def sensitivity = params.host_removal_verysensitive ? "--very-sensitive" : "--sensitive"
    def save_ids = params.host_removal_save_ids ? "Y" : "N"
    if (!meta.single_end){
        """
        bowtie2 -p ${task.cpus} \
                -x ${index[0].getSimpleName()} \
                -1 "${reads[0]}" -2 "${reads[1]}" \
                $sensitivity \
                --un-conc-gz ${prefix}.unmapped_%.fastq.gz \
                --al-conc-gz ${prefix}.mapped_%.fastq.gz \
                1> /dev/null \
                2> ${prefix}.bowtie2.log
        if [ ${save_ids} = "Y" ] ; then
            gunzip -c ${prefix}.mapped_1.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped_1.read_ids.txt
            gunzip -c ${prefix}.mapped_2.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped_2.read_ids.txt
        fi
        rm -f ${prefix}.mapped_*.fastq.gz

        echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
        """
    } else {
        """
        bowtie2 -p ${task.cpus} \
                -x ${index[0].getSimpleName()} \
                -U ${reads} \
                $sensitivity \
                --un-gz ${prefix}.unmapped.fastq.gz \
                --al-gz ${prefix}.mapped.fastq.gz \
                1> /dev/null \
                2> ${prefix}.bowtie2.log
        if [ ${save_ids} = "Y" ] ; then
            gunzip -c ${prefix}.mapped.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped.read_ids.txt
        fi
        rm -f ${prefix}.mapped.fastq.gz

        echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
        """
    }
}
