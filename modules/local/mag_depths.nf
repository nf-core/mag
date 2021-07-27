// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process MAG_DEPTHS {
    tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['assembler', 'id']) }

    // Using container from metabat2 process, since this will be anyway already downloaded and contains biopython and pandas
    conda (params.enable_conda ? "bioconda::metabat2=2.15 conda-forge::python=3.6.7 conda-forge::biopython=1.74 conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0"
    }

    input:
    tuple val(meta), path(bins)
    path(contig_depths)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.id}-binDepths.tsv"), emit: depths

    script:
    def software = getSoftwareName(task.process)
    """
    get_mag_depths.py --bins ${bins} \
                    --depths ${contig_depths} \
                    --assembly_name "${meta.assembler}-${meta.id}" \
                    --out "${meta.assembler}-${meta.id}-binDepths.tsv"
    """
}
