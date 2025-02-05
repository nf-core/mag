/*
 * Binning preparation with Bowtie2
 */

include { BOWTIE2_BUILD as BOWTIE2_ASSEMBLY_BUILD } from '../../modules/nf-core/bowtie2/build'
include { BOWTIE2_ALIGN as BOWTIE2_ASSEMBLY_ALIGN } from '../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_SORT                           } from '../../modules/nf-core/samtools/sort/main'

workflow BINNING_PREPARATION {
    take:
    assemblies // channel: [ val(meta), path(assembly) ]
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    // build bowtie2 index for all assemblies
    BOWTIE2_ASSEMBLY_BUILD(assemblies)

    // combine assemblies with sample reads for binning depending on specified mapping mode
    if (params.binning_map_mode == 'all') {
        // combine assemblies with reads of all samples
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.index.combine(reads)
    }
    else if (params.binning_map_mode == 'group') {
        // combine assemblies with reads of samples from same group
        ch_reads_bowtie2 = reads.map { meta, fastq -> [meta.group, meta, fastq] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.index
            .map { meta, assembly, index -> [meta.group, meta, assembly, index] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { _group, assembly_meta, assembly, index, reads_meta, fastq -> [assembly_meta, assembly, index, reads_meta, fastq] }
    }
    else {
        // i.e. --binning_map_mode 'own'
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map { meta, fastq -> [meta.id, meta, fastq] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.index
            .map { meta, assembly, index -> [meta.id, meta, assembly, index] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { _id, assembly_meta, assembly, index, reads_meta, fastq -> [assembly_meta, assembly, index, reads_meta, fastq] }
    }
    // TODO: run current mag with dumping `ch_bowtie2_input` before adn after alignment to see what info is in meta
    // and what is dropped. See if all should be exported, and either merge all meta into one one object or
    // call `meta2` etc. within `modules.config for prefix/tag

    ch_bowtie2_align_input = ch_bowtie2_input.map { id, assembly_meta, index, reads_meta, fastq ->
        [id, assembly_meta, index, reads_meta, fastq]
    }
    // TODO finish constructing this to go into BOWTIE2_ASSEMBLY_ALIGN

    BOWTIE2_ASSEMBLY_ALIGN(ch_bowtie2_input, [], [], false, false)
    // TODO fix ch_bowtie2_input to match, maybe single object?
    SAMTOOLS_SORT(BOWTIE2_ASSEMBLY_ALIGN.out.bam, [])
    ch_sorted_assembly_bams = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.csi)

    // group mappings for one assembly
    ch_grouped_mappings = ch_sorted_assembly_bams
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [meta, assembly.sort()[0], bams, bais] }

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { _assembly_meta, _reads_meta, log -> [log] }
    bowtie2_version          = BOWTIE2_ASSEMBLY_ALIGN.out.versions
    grouped_mappings         = ch_grouped_mappings
}
