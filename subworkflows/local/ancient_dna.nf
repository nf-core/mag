include { BCFTOOLS_CONSENSUS }                                                           from '../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PRE ; BCFTOOLS_INDEX as BCFTOOLS_INDEX_POST } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_VIEW }                                                                from '../../modules/nf-core/bcftools/view/main'
include { FREEBAYES }                                                                    from '../../modules/nf-core/freebayes/main'
include { PYDAMAGE_ANALYZE }                                                             from '../../modules/nf-core/pydamage/analyze/main'
include { PYDAMAGE_FILTER }                                                              from '../../modules/nf-core/pydamage/filter/main'
include { SAMTOOLS_FAIDX as FAIDX}                                                       from '../../modules/nf-core/samtools/faidx/main'

workflow ANCIENT_DNA_ASSEMBLY_VALIDATION {
    take:
        input //channel: [val(meta), path(contigs), path(bam), path(bam_index)]
    main:
        ch_versions = Channel.empty()

        PYDAMAGE_ANALYZE(
            input.map {
                meta, contigs, bam, bai -> [
                    meta, bam[0], bai[0]
                ]
            }
        )

        PYDAMAGE_FILTER(PYDAMAGE_ANALYZE.out.csv)
        ch_versions = ch_versions.mix(PYDAMAGE_ANALYZE.out.versions.first())

        if ( params.skip_ancient_damagecorrection ) {
            ch_corrected_contigs = Channel.empty()
        }

        if ( !params.skip_ancient_damagecorrection ) {
            FAIDX(input.map { item -> [ item[0], item[1] ] }, [[],[]] )
            freebayes_input = input.join(FAIDX.out.fai)  // [val(meta), path(contigs), path(bam), path(bam_index), path(fai)]
                                .multiMap{
                                    meta, contigs, bam, bai, fai ->
                                        reads: [ meta, bam, bai, [], [], [] ]
                                        fasta: [ contigs ]
                                        fai: [ fai ]
                                }
            FREEBAYES ( freebayes_input.reads,
                        freebayes_input.fasta,
                        freebayes_input.fai,
                        [],
                        [],
                        [] )

            BCFTOOLS_INDEX_PRE(FREEBAYES.out.vcf)
            BCFTOOLS_VIEW(FREEBAYES.out.vcf.join(BCFTOOLS_INDEX_PRE.out.tbi), [], [], [])
            BCFTOOLS_INDEX_POST(BCFTOOLS_VIEW.out.vcf)
            BCFTOOLS_CONSENSUS(BCFTOOLS_VIEW.out.vcf
                                    .join(BCFTOOLS_INDEX_POST.out.tbi)
                                    .join(input.map { item -> [ item[0], item[1] ] }))

            ch_corrected_contigs = BCFTOOLS_CONSENSUS.out.fasta

            ch_versions = ch_versions.mix(FAIDX.out.versions.first())
            ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())
            ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())
        }



    emit:
        contigs_recalled          = ch_corrected_contigs // channel: [ val(meta), path(fasta) ]
        pydamage_results          = PYDAMAGE_ANALYZE.out.csv     // channel: [ val(meta), path(csv) ]
        pydamage_filtered_results = PYDAMAGE_FILTER.out.csv      // channel: [ val(meta), path(csv) ]
        versions                  = ch_versions                  // channel: [ versions.yml ]
}

