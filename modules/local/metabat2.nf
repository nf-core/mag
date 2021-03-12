// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process METABAT2 {
    tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.assembler) }

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
    path "${meta.assembler}-${assembly}-depth.txt.gz"
    path "MetaBAT2/discarded/*"
    path '*.version.txt'                              , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    OMP_NUM_THREADS=${task.cpus} jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    gzip -c depth.txt > "${meta.assembler}-${assembly}-depth.txt.gz"
    metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o "MetaBAT2/${meta.assembler}-${meta.id}" -m ${params.min_contig_size} --unbinned --seed ${params.metabat_rng_seed}

    # save unbinned contigs above thresholds into individual files, dump others in one file
    split_fasta.py "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.fa" ${params.min_length_unbinned_contigs} ${params.max_unbinned_contigs} ${params.min_contig_size}

    mkdir MetaBAT2/discarded
    mv "MetaBAT2/${meta.assembler}-${meta.id}.lowDepth.fa" MetaBAT2/discarded/
    mv "MetaBAT2/${meta.assembler}-${meta.id}.tooShort.fa" MetaBAT2/discarded/
    mv "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.pooled.fa" MetaBAT2/discarded/
    mv "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.remaining.fa" MetaBAT2/discarded/

    # mv splitted file so that it doesnt end up in following processes
    mv "MetaBAT2/${meta.assembler}-${meta.id}.unbinned.fa" "${meta.assembler}-${meta.id}.unbinned.fa"

    echo \$(metabat2 --help 2>&1) | sed "s/^.*version 2\\://; s/ (Bioconda.*//" > ${software}.version.txt
    """
}
