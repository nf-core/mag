include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/main.nf'   addParams( options: [:] )
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PRE ; BCFTOOLS_INDEX as BCFTOOLS_INDEX_POST } from '../../modules/nf-core/modules/bcftools/index/main.nf'   addParams( options: [:] )
include { BCFTOOLS_VIEW } from '../../modules/nf-core/modules/bcftools/view/main.nf'   addParams( options: [:] )
include { FREEBAYES_GERMLINE as FREEBAYES } from '../../modules/nf-core/modules/freebayes/germline/main.nf'   addParams( options: [args : "-p 1"] )
include { PYDAMAGE_ANALYZE } from '../../modules/nf-core/modules/pydamage/analyze/main.nf'   addParams( options: [:] )
include { PYDAMAGE_FILTER } from '../../modules/nf-core/modules/pydamage/filter/main.nf'   addParams( options: [:] )
include { SAMTOOLS_FAIDX as FAIDX} from from '../../modules/nf-core/modules/samtools/faidx/main.nf'   addParams( options: [:] )

worflow ANCIENT_DNA_ASSEMLY_VALIDATION {
    take:
        contigs // assembled contigs [ val(meta), [contigs] ]
        bam     // reads mapped back on contigs [ val(meta), [bam] ]
        bam_index // bam index [ val(meta), [index] ]
    main:
        PYDAMAGE_ANALYZE(contigs
                            .join(bam)
                            .join(bam_index))
        PYDAMAGE_FILTER(PYDAMAGE_ANALYZE.out.csv)
        FAIDX(contigs.map { it -> [ it[1] ] })
        FREEBAYES (bam.join(bam_index), contigs.map { it -> [ it[1] ] }, FAIDX.out.FAI )
        BCFTOOLS_INDEX_PRE(FREEBAYES.out.vcf)
        BCFTOOLS_VIEW(FREEBAYES.out.vcf.join(BCFTOOLS_INDEX_PRE.out.tbi)) // Add filtering parameters
        BCFTOOLS_INDEX_POST(BCFTOOLS_VIEW.out.vcf)
        BCFTOOLS_CONSENSUS(BCFTOOLS_VIEW.out.vcf
                               .join(BCFTOOLS_INDEX_POST.out.tbi)
                               .join(contigs))
    emit:
        contigs_recalled: BCFTOOLS_CONSENSUS.out.fasta
        pydamage_results: PYDAMAGE_ANALYZE.out.csv
        pydamage_filtered_results: PYDAMAGE_FILTER.out.csv
        


}

