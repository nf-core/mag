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
    ch_versions = Channel.empty()

    // build bowtie2 index for all assemblies
    BOWTIE2_ASSEMBLY_BUILD(assemblies)
    ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_BUILD.out.versions.first())


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
            .map { _group, assembly_meta, assembly, index, reads_meta, fastq ->
                [assembly_meta, assembly, index, reads_meta, fastq]
            }
    }
    else {
        // i.e. --binning_map_mode 'own'
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map { meta, fastq -> [meta.id, meta, fastq] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.index
            .map { meta, assembly, index -> [meta.id, meta, assembly, index] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { _id, assembly_meta, assembly, index, reads_meta, fastq ->
                [assembly_meta, assembly, index, reads_meta, fastq]
            }
    }

    // RECONFIGURE CHANNEL TO MATCH BOWTIE2_ALIGN INPUT STRUCTURE
    // I purposely flip the location of the metas here, to ensure the assembly meta gets
    // exported with the resulting BAM file from the bowtie2 module, and not the reads meta.
    // We have to retain the reads meta for use within prefix and tag in the bowtie2 module.
    // I don't combine assembly_meta and relevant reads_meta, as otherwise would have a map
    // operation after each use of SAMTOOLS_SORT.out.*, to remove the now irrelevant reads_meta
    // information.
    ch_bowtie2_align_input = ch_bowtie2_input.multiMap { assembly_meta, assembly, index, reads_meta, fastq ->
        reads: [assembly_meta, fastq]
        index: [[:], index]
        assembly: [reads_meta, assembly]
        saveunaligned: false
        sortbam: false
    }

    BOWTIE2_ASSEMBLY_ALIGN(
        ch_bowtie2_align_input.reads,
        ch_bowtie2_align_input.index,
        ch_bowtie2_align_input.ref,
        ch_bowtie2_align_input.saveunaligned,
        ch_bowtie2_align_input.sortbam,
    )
    ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_ALIGN.out.versions.first())

    SAMTOOLS_SORT(BOWTIE2_ASSEMBLY_ALIGN.out.bam, [[:], []])
    ch_sorted_assembly_bams = assemblies
        .join(SAMTOOLS_SORT.out.bam)
        .join(SAMTOOLS_SORT.out.csi)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // group mappings for one assembly
    ch_grouped_mappings = ch_sorted_assembly_bams
        .groupTuple(by: 0)
        .map { meta, assembly, bams, csis -> [meta, assembly.sort()[0], bams, csis] }
        .dump(tag: 'ch_grouped_mappings')

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { _meta, log -> [log] }
    version                  = ch_versions
    grouped_mappings         = ch_grouped_mappings
}
