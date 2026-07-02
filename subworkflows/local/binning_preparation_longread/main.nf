include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ASSEMBLY_ALIGN } from '../../../modules/nf-core/minimap2/align/main'

def mappingKey(meta) {
    if (params.binning_map_mode == 'group') {
        return meta.group
    }
    if (params.binning_map_mode == 'own') {
        return meta.id
    }
    return 'all'
}

workflow LONGREAD_BINNING_PREPARATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_reads // [val(meta), path(reads)]

    main:
    ch_versions = channel.empty()

    MINIMAP2_ASSEMBLY_INDEX(ch_assemblies)

    // combine assemblies with reads depending on mapping mode
    ch_reads_minimap2 = ch_reads.map { meta, reads -> [mappingKey(meta), meta, reads] }
    ch_minimap2_input = MINIMAP2_ASSEMBLY_INDEX.out.index
        .map { meta_idx, index -> [mappingKey(meta_idx), meta_idx, index] }
        .combine(ch_reads_minimap2, by: 0)
        .multiMap { _key, meta_idx, index, meta_reads, reads ->
            reads: [meta_idx, reads]
            index: [meta_reads, index]
        }

    MINIMAP2_ASSEMBLY_ALIGN(ch_minimap2_input.reads, ch_minimap2_input.index, true, 'bai', false, false)

    // read sets expected per assembly, counted from the reads so groupTuple can
    // release each assembly as soon as its own alignments finish instead of
    // waiting for every MINIMAP2_ASSEMBLY_ALIGN task
    ch_reads_count = ch_reads
        .map { meta, _reads -> [mappingKey(meta), 1] }
        .groupTuple()
        .map { key, counts -> [key, counts.size()] }

    ch_mapping_counts = ch_assemblies
        .map { meta, _assembly -> [mappingKey(meta), meta.assembler, meta.id] }
        .combine(ch_reads_count, by: 0)
        .map { _key, assembler, id, count -> [[assembler, id], count] }

    ch_grouped_mappings = MINIMAP2_ASSEMBLY_ALIGN.out.bam
        .join(MINIMAP2_ASSEMBLY_ALIGN.out.index)
        .map { meta, bam, bai -> [[meta.assembler, meta.id], meta, bam, bai] }
        .combine(ch_mapping_counts, by: 0)
        .map { key, meta, bam, bai, count -> [groupKey(key, count), meta, bam, bai] }
        .groupTuple()
        .map { _key, metas, bams, bais -> [metas[0], bams, bais] }
        .combine(ch_assemblies, by: 0)
        .map { meta, bams, bais, assembly -> [meta, assembly, bams, bais] }

    emit:
    versions         = ch_versions
    grouped_mappings = ch_grouped_mappings
}
