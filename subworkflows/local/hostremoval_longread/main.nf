//
// Remove host reads via alignment and export off-target reads
//

include { MINIMAP2_INDEX as MINIMAP2_HOST_INDEX              } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_HOST_ALIGN              } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_HOSTREMOVED_INDEX       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as SAMTOOLS_HOSTREMOVED_STATS       } from '../../../modules/nf-core/samtools/stats/main'

include { SAMTOOLS_UNMAPPED as SAMTOOLS_HOSTREMOVED_UNMAPPED } from '../../../modules/local/samtools_unmapped/main'


workflow LONGREAD_HOSTREMOVAL {
    take:
    ch_reads      // [val(meta), path(fastq)]
    val_reference // path

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_host_reference = val_reference.map { ref -> [[:], ref] }
    ch_host_fasta_for_build = ch_host_reference
        .combine(ch_reads)
        .map { host_meta, host_fasta, _meta, _reads ->
            [host_meta, host_fasta]
        }
        .first()
    // makes sure to only use the host fasta if the long read channel is not empty
    ch_minimap2_index = MINIMAP2_HOST_INDEX(ch_host_fasta_for_build).index
    ch_versions = ch_versions.mix(MINIMAP2_HOST_INDEX.out.versions)

    MINIMAP2_HOST_ALIGN(ch_reads, ch_minimap2_index, true, 'bai', false, false)
    ch_minimap2_mapped = MINIMAP2_HOST_ALIGN.out.bam.map { meta, reads ->
        [meta, reads, []]
    }
    ch_versions = ch_versions.mix(MINIMAP2_HOST_ALIGN.out.versions)

    // Generate unmapped reads FASTQ for downstream taxprofiling
    SAMTOOLS_HOSTREMOVED_UNMAPPED(ch_minimap2_mapped)
    ch_versions = ch_versions.mix(SAMTOOLS_HOSTREMOVED_UNMAPPED.out.versions)

    // Indexing whole BAM for host removal statistics
    SAMTOOLS_HOSTREMOVED_INDEX(MINIMAP2_HOST_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_HOSTREMOVED_INDEX.out.versions)

    bam_bai = MINIMAP2_HOST_ALIGN.out.bam.join(SAMTOOLS_HOSTREMOVED_INDEX.out.bai)

    SAMTOOLS_HOSTREMOVED_STATS(bam_bai, ch_host_fasta_for_build)
    ch_versions = ch_versions.mix(SAMTOOLS_HOSTREMOVED_STATS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_HOSTREMOVED_STATS.out.stats)

    emit:
    stats         = SAMTOOLS_HOSTREMOVED_STATS.out.stats //channel: [val(meta), [reads  ] ]
    reads         = SAMTOOLS_HOSTREMOVED_UNMAPPED.out.fastq // channel: [ val(meta), [ reads ] ]
    versions      = ch_versions // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
