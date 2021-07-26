// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process METABAT2 {
    tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['assembler']) }

    conda (params.enable_conda ? "bioconda::metabat2=2.15 conda-forge::python=3.6.7 conda-forge::biopython=1.74 conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0"
    }

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("MetaBAT2/*.fa")            , emit: bins
    path "${meta.assembler}-${meta.id}-depth.txt.gz"  , emit: depths
    path "MetaBAT2/discarded/*"                       , emit: discarded
    path '*.version.txt'                              , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    OMP_NUM_THREADS=${task.cpus} jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o "MetaBAT2/${meta.assembler}-${meta.id}" -m ${params.min_contig_size} --unbinned --seed ${params.metabat_rng_seed}

    gzip depth.txt
    mv depth.txt.gz "${meta.assembler}-${meta.id}-depth.txt.gz"

    # save unbinned contigs above thresholds into individual files, dump others in one file
    split_fasta.py "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.fa" ${params.min_length_unbinned_contigs} ${params.max_unbinned_contigs} ${params.min_contig_size}

    # delete splitted file so that it doesnt end up in following processes
    rm "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.fa"

    mkdir MetaBAT2/discarded
    gzip "MetaBAT2/${meta.assembler}-${meta.id}.lowDepth.fa" \
        "MetaBAT2/${meta.assembler}-${meta.id}.tooShort.fa" \
        "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.pooled.fa" \
        "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.remaining.fa"
    mv "MetaBAT2/${meta.assembler}-${meta.id}".*.fa.gz MetaBAT2/discarded/

    echo \$(metabat2 --help 2>&1) | sed "s/^.*version 2\\://; s/ (Bioconda.*//" > ${software}.version.txt
    """
}
