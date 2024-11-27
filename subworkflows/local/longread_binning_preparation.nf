include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX                 } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ASSEMBLY_ALIGN                 } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_ASSEMBLY_INDEX                 } from '../../modules/nf-core/samtools/index/main'

workflow LONGREAD_BINNING_PREPARATION {
    take:
    assemblies // channel: [ val(meta), path(assembly) ]
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    MINIMAP2_ASSEMBLY_INDEX ( assemblies )
    ch_versions       = ch_versions.mix( MINIMAP2_ASSEMBLY_INDEX.out.versions )

    // combine assemblies with sample reads for binning depending on specified mapping mode
    if (params.binning_map_mode == 'all'){
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .combine(reads)
    } else if (params.binning_map_mode == 'group'){
        ch_reads_minimap2 = reads.map{ meta, reads -> [ meta.group, meta, reads ] }
        ch_index_minimap2 = MINIMAP2_ASSEMBLY_INDEX.out.index
            .combine(assemblies)
            .map {index, assembly -> [assembly[0], assembly[1], index[1]] }
        ch_minimap2_input = ch_index_minimap2
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .combine(ch_reads_minimap2, by: 0)
            .map { group, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }

    } else {
        // i.e. --binning_map_mode 'own'
        ch_reads_minimap2 = reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_index_minimap2 = MINIMAP2_ASSEMBLY_INDEX.out.index
            .combine(assemblies)
            .map {index, assembly -> [assembly[0], assembly[1], index[1]] }
        ch_minimap2_input = ch_index_minimap2
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_minimap2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    }

    MINIMAP2_ASSEMBLY_ALIGN ( ch_minimap2_input )
    SAMTOOLS_ASSEMBLY_INDEX ( MINIMAP2_ASSEMBLY_ALIGN.out.bam )
    ch_versions      = ch_versions.mix( MINIMAP2_ASSEMBLY_ALIGN.out.versions.first() )
    ch_versions      = ch_versions.mix( SAMTOOLS_ASSEMBLY_INDEX.out.versions.first() )

    ch_minimap2_aligned = MINIMAP2_ASSEMBLY_ALIGN.out.bam
        .combine(SAMTOOLS_ASSEMBLY_INDEX.out.bai, by: 0)

    ch_grouped_mappings = ch_minimap2_input
        .map{assembly_meta, assembly, index, reads_meta, reads -> [assembly_meta, assembly]}
        .combine(ch_minimap2_aligned, by: 0)
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }

    emit:
    versions           = ch_versions
    grouped_mappings   = ch_grouped_mappings




    // group mappings for one assembly


}
