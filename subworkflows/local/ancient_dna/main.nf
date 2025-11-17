include { BCFTOOLS_CONSENSUS                    } from '../../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_INDEX                        } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_VIEW                         } from '../../../modules/nf-core/bcftools/view/main'
include { FREEBAYES                             } from '../../../modules/nf-core/freebayes/main'
include { PYDAMAGE_ANALYZE                      } from '../../../modules/nf-core/pydamage/analyze/main'
include { PYDAMAGE_FILTER                       } from '../../../modules/nf-core/pydamage/filter/main'
include { SAMTOOLS_FAIDX as FAIDX               } from '../../../modules/nf-core/samtools/faidx/main'

workflow ANCIENT_DNA_ASSEMBLY_VALIDATION {
    take:
    ch_input // [val(meta), path(contigs), path(bam), path(bam_index)]

    main:
    ch_versions = channel.empty()

    PYDAMAGE_ANALYZE(
        ch_input.map { meta, _contigs, bam, bai ->
            [
                meta,
                bam[0],
                bai[0],
            ]
        }
    )
    ch_versions = ch_versions.mix(PYDAMAGE_ANALYZE.out.versions)

    PYDAMAGE_FILTER(PYDAMAGE_ANALYZE.out.csv)
    ch_versions = ch_versions.mix(PYDAMAGE_FILTER.out.versions)

    if (params.skip_ancient_damagecorrection) {
        ch_corrected_contigs = channel.empty()
    }

    if (!params.skip_ancient_damagecorrection) {
        FAIDX(ch_input.map { item -> [item[0], item[1]] }, [[], []], false)
        ch_versions = ch_versions.mix(FAIDX.out.versions)

        freebayes_input = ch_input
            .join(FAIDX.out.fai)
            .multiMap { meta, contigs, bam, bai, fai ->
                reads: [meta, bam, bai, [], [], []]
                fasta: [meta, contigs]
                fai: [meta, fai]
            }
        FREEBAYES(
            freebayes_input.reads,
            freebayes_input.fasta,
            freebayes_input.fai,
            [[], []],
            [[], []],
            [[], []],
        )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)

        BCFTOOLS_INDEX(FREEBAYES.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
        BCFTOOLS_VIEW(FREEBAYES.out.vcf.join(BCFTOOLS_INDEX.out.tbi), [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        BCFTOOLS_CONSENSUS(
            BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi).join(ch_input.map { item -> [item[0], item[1], []] })
        )
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

        ch_corrected_contigs = BCFTOOLS_CONSENSUS.out.fasta
    }

    emit:
    contigs_recalled          = ch_corrected_contigs // channel: [ val(meta), path(fasta) ]
    pydamage_results          = PYDAMAGE_ANALYZE.out.csv // channel: [ val(meta), path(csv) ]
    pydamage_filtered_results = PYDAMAGE_FILTER.out.csv // channel: [ val(meta), path(csv) ]
    versions                  = ch_versions // channel: [ versions.yml ]
}
