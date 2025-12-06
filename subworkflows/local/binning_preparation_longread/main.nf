include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ASSEMBLY_ALIGN } from '../../../modules/nf-core/minimap2/align/main'

workflow LONGREAD_BINNING_PREPARATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_reads      // [val(meta), path(reads)]

    main:
    ch_versions = channel.empty()

    MINIMAP2_ASSEMBLY_INDEX(ch_assemblies)
    ch_versions = ch_versions.mix(MINIMAP2_ASSEMBLY_INDEX.out.versions)

    if (params.binning_map_mode == 'all') {
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .combine(ch_reads)
            .map { meta_idx, idx, meta_reads, reads -> [meta_idx, idx, meta_reads, reads] }
    }
    else if (params.binning_map_mode == 'group') {
        ch_reads_minimap2 = ch_reads.map { meta, reads -> [meta.group, meta, reads] }
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .map { meta_idx, index -> [meta_idx.group, meta_idx, index] }
            .combine(ch_reads_minimap2, by: 0)
            .map { _group, meta_idx, idx, meta_reads, reads -> [meta_idx, idx, meta_reads, reads] }
    }
    else {
        ch_reads_minimap2 = ch_reads.map { meta, reads -> [meta.id, meta, reads] }
        ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
            .map { meta_idx, index -> [meta_idx.id, meta_idx, index] }
            .combine(ch_reads_minimap2, by: 0)
            .map { _id, meta_idx, idx, meta_reads, reads -> [meta_idx, idx, meta_reads, reads] }
    }

    ch_minimap2_input_reads = ch_minimap2_input.map { meta_idx, _index, _meta, reads -> [meta_idx, reads] }
    ch_minimap2_input_idx = ch_minimap2_input.map { _meta_idx, index, meta, _reads -> [meta, index] }

    MINIMAP2_ASSEMBLY_ALIGN(ch_minimap2_input_reads, ch_minimap2_input_idx, true, 'bai', false, false)
    ch_versions = ch_versions.mix(MINIMAP2_ASSEMBLY_ALIGN.out.versions)

    ch_grouped_mappings_reads = MINIMAP2_ASSEMBLY_ALIGN.out.bam.groupTuple(by: 0)
    ch_grouped_mappings_index = MINIMAP2_ASSEMBLY_ALIGN.out.index.groupTuple(by: 0)
    ch_grouped_mappings = ch_grouped_mappings_reads
        .combine(ch_grouped_mappings_index, by: 0)
        .combine(ch_assemblies, by: 0)
        .map { meta, bams, bais, assembly -> [meta, assembly, bams, bais] }

    emit:
    versions         = ch_versions
    grouped_mappings = ch_grouped_mappings
}
