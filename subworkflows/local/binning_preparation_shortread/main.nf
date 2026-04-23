/*
 * Binning preparation for short reads
 */

include { BOWTIE2_ASSEMBLY_BUILD } from '../../../modules/local/bowtie2_assembly_build/main'
include { BOWTIE2_ASSEMBLY_ALIGN } from '../../../modules/local/bowtie2_assembly_align/main'
include { COVERM_CONTIG as COVERM_CONTIG_SHORTREAD } from '../../../modules/nf-core/coverm/contig/main'
include { SAMTOOLS_INDEX as SAMTOOLS_COVERM_SHORTREAD_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow SHORTREAD_BINNING_PREPARATION {
    take:
    ch_assemblies // [val(meta), path(assembly)]
    ch_reads      // [val(meta), path(reads)]

    main:

    ch_versions = channel.empty()
    ch_multiqc  = channel.empty()
    ch_contig_depths = channel.empty()

    if (params.shortread_coverage_mapper == 'bowtie2') {
        // build bowtie2 index for all assemblies
        BOWTIE2_ASSEMBLY_BUILD(ch_assemblies)
        ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_BUILD.out.versions)

        // combine assemblies with sample reads for binning depending on specified mapping mode
        if (params.binning_map_mode == 'all') {
            // combine assemblies with reads of all samples
            ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index.combine(ch_reads)
        }
        else if (params.binning_map_mode == 'group') {
            // combine assemblies with reads of samples from same group
            ch_reads_bowtie2 = ch_reads.map { meta, sample_reads -> [meta.group, meta, sample_reads] }
            ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
                .map { meta, assembly, index -> [meta.group, meta, assembly, index] }
                .combine(ch_reads_bowtie2, by: 0)
                .map { _group, assembly_meta, assembly, index, reads_meta, sample_reads ->
                    [assembly_meta, assembly, index, reads_meta, sample_reads]
                }
        }
        else {
            // i.e. --binning_map_mode 'own'
            // combine assemblies (not co-assembled) with reads from own sample
            ch_reads_bowtie2 = ch_reads.map { meta, sample_reads -> [meta.id, meta, sample_reads] }
            ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
                .map { meta, assembly, index -> [meta.id, meta, assembly, index] }
                .combine(ch_reads_bowtie2, by: 0)
                .map { _id, assembly_meta, assembly, index, reads_meta, sample_reads ->
                    [assembly_meta, assembly, index, reads_meta, sample_reads]
                }
        }

        BOWTIE2_ASSEMBLY_ALIGN(ch_bowtie2_input)
        ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_ALIGN.out.versions)

        // group mappings for one assembly
        ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGN.out.mappings
            .groupTuple(by: 0)
            .map { meta, assembly, bams, bais -> [meta, assembly.sort()[0], bams, bais] }

        ch_coverm_input = ch_grouped_mappings.multiMap { meta, _assembly, bams, _bais ->
            reads:     [meta, bams]
            reference: [meta, []]
        }
        COVERM_CONTIG_SHORTREAD(ch_coverm_input.reads, ch_coverm_input.reference, true, false, false)
        ch_contig_depths = COVERM_CONTIG_SHORTREAD.out.coverage

        ch_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { _assembly_meta, _reads_meta, log -> [log] }
    }
    else {
        // CoverM native mapper: skip bowtie2 alignment, group reads with assemblies
        // so CoverM can handle alignment and coverage calculation internally
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
            // i.e. --binning_map_mode 'own'
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
            reads:     [meta, reads]
            reference: [meta, assembly]
        }
        COVERM_CONTIG_SHORTREAD(ch_coverm_input_multi.reads, ch_coverm_input_multi.reference, false, false, true)
        ch_contig_depths = COVERM_CONTIG_SHORTREAD.out.coverage

        SAMTOOLS_COVERM_SHORTREAD_INDEX(COVERM_CONTIG_SHORTREAD.out.bam.transpose())
        ch_versions = ch_versions.mix(SAMTOOLS_COVERM_SHORTREAD_INDEX.out.versions)

        ch_grouped_mappings = COVERM_CONTIG_SHORTREAD.out.bam
            .combine(SAMTOOLS_COVERM_SHORTREAD_INDEX.out.bai.groupTuple(by: 0), by: 0)
            .combine(ch_coverm_input.map { meta, assembly, _reads -> [meta, assembly] }, by: 0)
            .map { meta, bams, bais, assembly -> [meta, assembly, bams, bais] }
    }

    emit:
    bowtie2_assembly_multiqc = ch_multiqc
    versions                 = ch_versions
    grouped_mappings         = ch_grouped_mappings
    contig_depths            = ch_contig_depths
}
