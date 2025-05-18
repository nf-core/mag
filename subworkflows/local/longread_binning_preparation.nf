include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX                 } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ASSEMBLY_ALIGN                 } from '../../modules/nf-core/minimap2/align/main'

workflow LONGREAD_BINNING_PREPARATION {
    take:
    assemblies // channel: [ val(meta), path(assembly) ]
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions       = Channel.empty()

    MINIMAP2_ASSEMBLY_INDEX ( assemblies )

    if (params.binning_map_mode == 'all'){
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .combine(reads)
            .map { meta_idx, idx, meta_reads, reads -> [ meta_idx, idx,  meta_reads, reads ] }

    } else if (params.binning_map_mode == 'group') {
        ch_reads_minimap2 = reads.map{ meta, reads -> [ meta.group, meta, reads ] }
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .map { meta_idx, index -> [ meta_idx.group, meta_idx, index ] }
            .combine(ch_reads_minimap2, by: 0)
            .map { group, meta_idx, idx, meta_reads, reads -> [ meta_idx, idx, meta_reads, reads ] }

    } else {
        ch_reads_minimap2 = reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .map { meta_idx, index -> [ meta_idx.id, meta_idx, index ] }
            .combine(ch_reads_minimap2, by: 0)
            .map { id, meta_idx, idx, meta_reads, reads -> [ meta_idx, idx, meta_reads, reads ] }
    }

    ch_minimap2_input_reads = ch_minimap2_input
        .map { meta_idx, index, meta, reads -> [ meta_idx, reads ] }
    ch_minimap2_input_idx = ch_minimap2_input
        .map { meta_idx, index, meta, reads -> [ meta, index ] }

    MINIMAP2_ASSEMBLY_ALIGN ( ch_minimap2_input_reads, ch_minimap2_input_idx, true, 'bai', false, false )
    ch_versions      = ch_versions.mix( MINIMAP2_ASSEMBLY_ALIGN.out.versions.first() )

    ch_grouped_mappings_reads = MINIMAP2_ASSEMBLY_ALIGN.out.bam
        .groupTuple(by: 0)
    ch_grouped_mappings_index = MINIMAP2_ASSEMBLY_ALIGN.out.index
        .groupTuple(by: 0)
    ch_grouped_mappings = ch_grouped_mappings_reads
        .combine(ch_grouped_mappings_index, by: 0)
        .combine(assemblies, by: 0)
        .map { meta, bams, bais, assembly -> [ meta, assembly, bams, bais ] }

    emit:
    versions           = ch_versions
    grouped_mappings   = ch_grouped_mappings
}
