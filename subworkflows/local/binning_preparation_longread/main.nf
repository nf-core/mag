include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ASSEMBLY_ALIGN } from '../../../modules/nf-core/minimap2/align/main'

workflow LONGREAD_BINNING_PREPARATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_reads // [val(meta), path(reads)]

    main:
    ch_versions = channel.empty()

    MINIMAP2_ASSEMBLY_INDEX(ch_assemblies)

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

    // expected number of mappings (read sets) per assembly, resolved early from
    // the alignment inputs so groupTuple can release each assembly's group as
    // soon as it is complete instead of waiting for every MINIMAP2_ASSEMBLY_ALIGN
    // task to finish
    ch_mapping_counts = ch_minimap2_input
        .map { meta_idx, _index, _meta_reads, _reads -> [[meta_idx.assembler, meta_idx.id], 1] }
        .groupTuple()
        .map { key, counts -> [key, counts.size()] }

    ch_grouped_mappings_reads = MINIMAP2_ASSEMBLY_ALIGN.out.bam
        .map { meta, bam -> [[meta.assembler, meta.id], meta, bam] }
        .combine(ch_mapping_counts, by: 0)
        .map { key, meta, bam, count -> [groupKey(key, count), meta, bam] }
        .groupTuple()
        .map { _key, metas, bams -> [metas[0], bams] }
    ch_grouped_mappings_index = MINIMAP2_ASSEMBLY_ALIGN.out.index
        .map { meta, bai -> [[meta.assembler, meta.id], meta, bai] }
        .combine(ch_mapping_counts, by: 0)
        .map { key, meta, bai, count -> [groupKey(key, count), meta, bai] }
        .groupTuple()
        .map { _key, metas, bais -> [metas[0], bais] }
    ch_grouped_mappings = ch_grouped_mappings_reads
        .combine(ch_grouped_mappings_index, by: 0)
        .combine(ch_assemblies, by: 0)
        .map { meta, bams, bais, assembly -> [meta, assembly, bams, bais] }

    emit:
    versions         = ch_versions
    grouped_mappings = ch_grouped_mappings
}
