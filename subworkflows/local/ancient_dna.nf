include { BCFTOOLS_CONSENSUS }                                                           from '../../modules/nf-core/modules/bcftools/consensus/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PRE ; BCFTOOLS_INDEX as BCFTOOLS_INDEX_POST } from '../../modules/nf-core/modules/bcftools/index/main'
include { BCFTOOLS_VIEW }                                                                from '../../modules/nf-core/modules/bcftools/view/main'
include { FREEBAYES }                                                                    from '../../modules/nf-core/modules/freebayes/main'
include { PYDAMAGE_ANALYZE }                                                             from '../../modules/nf-core/modules/pydamage/analyze/main'
include { PYDAMAGE_FILTER }                                                              from '../../modules/nf-core/modules/pydamage/filter/main'
include { SAMTOOLS_FAIDX as FAIDX}                                                       from '../../modules/nf-core/modules/samtools/faidx/main'

workflow ANCIENT_DNA_ASSEMLY_VALIDATION {
    take:
        input //channel: [val(meta), path(contigs), path(bam), path(bam_index)]
    main:
        PYDAMAGE_ANALYZE(input.map {item -> [item[0], item[2], item[3]]})
        PYDAMAGE_FILTER(PYDAMAGE_ANALYZE.out.csv)
        FAIDX(input.map { item -> [ item[0], item[1] ] })
        input.join(FAIDX.out.fai)
            .set { freebayes_input } // [val(meta), path(contigs), path(bam), path(bam_index), path(fai)]
        FREEBAYES (freebayes_input.map { item -> [item[0], item[2], item[3], [], [], []] },
                    freebayes_input.map { item -> item[1] },
                    freebayes_input.map { item -> item[4] },
                    [],
                    [],
                    [] )

        BCFTOOLS_INDEX_PRE(FREEBAYES.out.vcf)
        BCFTOOLS_VIEW(FREEBAYES.out.vcf.join(BCFTOOLS_INDEX_PRE.out.tbi), [], [], [])
        BCFTOOLS_INDEX_POST(BCFTOOLS_VIEW.out.vcf)
        BCFTOOLS_CONSENSUS(BCFTOOLS_VIEW.out.vcf
                                .join(BCFTOOLS_INDEX_POST.out.tbi)
                                .join(input.map { item -> [ item[0], item[1] ] }))

        ch_versions = Channel.empty()
        ch_versions = PYDAMAGE_ANALYZE.out.versions.first()
        ch_versions = ch_versions.mix(FAIDX.out.versions.first())
        ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())
    emit:
        contigs_recalled          = BCFTOOLS_CONSENSUS.out.fasta // channel: [ val(meta), path(fasta) ]
        pydamage_results          = PYDAMAGE_ANALYZE.out.csv     // channel: [ val(meta), path(csv) ]
        pydamage_filtered_results = PYDAMAGE_FILTER.out.csv      // channel: [ val(meta), path(csv) ]
        versions                  = ch_versions                  // channel: [ versions.yml ]
}

