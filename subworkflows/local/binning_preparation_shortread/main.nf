/*
 * Binning preparation with Bowtie2
 */

include { BOWTIE2_ASSEMBLY_BUILD } from '../../../modules/local/bowtie2_assembly_build/main'
include { BOWTIE2_ASSEMBLY_ALIGN } from '../../../modules/local/bowtie2_assembly_align/main'

def mappingKey(meta) {
    if (params.binning_map_mode == 'group') {
        return meta.group
    }
    if (params.binning_map_mode == 'own') {
        return meta.id
    }
    return 'all'
}

workflow SHORTREAD_BINNING_PREPARATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_reads      // [val(meta), path(reads)]

    main:

    ch_versions = channel.empty()
    // build bowtie2 index for all assemblies
    BOWTIE2_ASSEMBLY_BUILD(ch_assemblies)
    ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_BUILD.out.versions)

    // combine assemblies with sample reads depending on mapping mode
    ch_reads_bowtie2 = ch_reads.map { meta, sample_reads -> [mappingKey(meta), meta, sample_reads] }
    ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
        .map { meta, assembly, index -> [mappingKey(meta), meta, assembly, index] }
        .combine(ch_reads_bowtie2, by: 0)
        .map { _key, assembly_meta, assembly, index, reads_meta, sample_reads ->
            [assembly_meta, assembly, index, reads_meta, sample_reads]
        }

    BOWTIE2_ASSEMBLY_ALIGN(ch_bowtie2_input)
    ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_ALIGN.out.versions)

    // read sets expected per assembly, counted from the reads so groupTuple can
    // release each assembly as soon as its own alignments finish instead of
    // waiting for every BOWTIE2_ASSEMBLY_ALIGN task
    ch_reads_count = ch_reads
        .map { meta, _reads -> [mappingKey(meta), 1] }
        .groupTuple()
        .map { key, counts -> [key, counts.size()] }

    ch_mapping_counts = ch_assemblies
        .map { meta, _assembly -> [mappingKey(meta), meta.assembler, meta.id] }
        .combine(ch_reads_count, by: 0)
        .map { _key, assembler, id, count -> [[assembler, id], count] }

    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGN.out.mappings
        .map { meta, assembly, bam, bai -> [[meta.assembler, meta.id], meta, assembly, bam, bai] }
        .combine(ch_mapping_counts, by: 0)
        .map { key, meta, assembly, bam, bai, count -> [groupKey(key, count), meta, assembly, bam, bai] }
        .groupTuple()
        .map { _key, metas, assemblies, bams, bais -> [metas[0], assemblies.toSorted()[0], bams, bais] }

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { _assembly_meta, _reads_meta, log -> [log] }
    versions                 = ch_versions
    grouped_mappings         = ch_grouped_mappings
}
