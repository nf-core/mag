include { COVERM_CONTIG as COVERM_CONTIG_LONGREAD } from '../../../modules/nf-core/coverm/contig/main'
include { SAMTOOLS_INDEX as SAMTOOLS_COVERM_LONGREAD_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow LONGREAD_BINNING_PREPARATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_reads // [val(meta), path(reads)]

    main:
    ch_versions = channel.empty()
    ch_contig_depths = channel.empty()

    if (params.binning_map_mode == 'all') {
        ch_coverm_reads = ch_assemblies
            .combine(ch_reads)
            .map { assembly_meta, assembly, _reads_meta, reads -> [assembly_meta, assembly, reads] }
    }
    else if (params.binning_map_mode == 'group') {
        ch_reads_grouped = ch_reads.map { meta, reads -> [meta.group, reads] }
        ch_coverm_reads  = ch_assemblies
            .map { meta, assembly -> [meta.group, meta, assembly] }
            .combine(ch_reads_grouped, by: 0)
            .map { _group, assembly_meta, assembly, reads -> [assembly_meta, assembly, reads] }
    }
    else {
        ch_reads_grouped = ch_reads.map { meta, reads -> [meta.id, reads] }
        ch_coverm_reads  = ch_assemblies
            .map { meta, assembly -> [meta.id, meta, assembly] }
            .combine(ch_reads_grouped, by: 0)
            .map { _id, assembly_meta, assembly, reads -> [assembly_meta, assembly, reads] }
    }

    ch_coverm_input = ch_coverm_reads
        .groupTuple(by: 0)
        .map { meta, assembly, reads_list ->
            [meta, assembly.sort()[0], reads_list.flatten()]
        }

    ch_coverm_input_multi = ch_coverm_input.multiMap { meta, assembly, reads ->
        reads:     [meta + [single_end: true], reads]
        reference: [meta + [single_end: true], assembly]
    }
    COVERM_CONTIG_LONGREAD(ch_coverm_input_multi.reads, ch_coverm_input_multi.reference, false, false, true)
    ch_contig_depths = COVERM_CONTIG_LONGREAD.out.coverage.map { meta, depth -> [meta - meta.subMap('single_end'), depth] }
    ch_coverm_bam = COVERM_CONTIG_LONGREAD.out.bam.map { meta, bams -> [meta - meta.subMap('single_end'), bams] }

    SAMTOOLS_COVERM_LONGREAD_INDEX(ch_coverm_bam.transpose())
    ch_versions = ch_versions.mix(SAMTOOLS_COVERM_LONGREAD_INDEX.out.versions)

    ch_grouped_mappings = ch_coverm_bam
        .combine(SAMTOOLS_COVERM_LONGREAD_INDEX.out.bai.groupTuple(by: 0), by: 0)
        .combine(ch_coverm_input.map { meta, assembly, _reads -> [meta, assembly] }, by: 0)
        .map { meta, bams, bais, assembly -> [meta, assembly, bams, bais] }

    emit:
    versions         = ch_versions
    grouped_mappings = ch_grouped_mappings
    contig_depths    = ch_contig_depths
}
