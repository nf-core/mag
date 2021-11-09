params.bcftools_consensus_options = [:]
params.bcftools_view_options      = [:]
params.bcftools_index_options     = [:]
params.freebayes_options          = [:]
params.pydamage_analyze_options   = [:]
params.pydamage_filter_options    = [:]
params.samtools_faidx_options     = [:]

include { BCFTOOLS_CONSENSUS }                                                           from '../../modules/nf-core/modules/bcftools/consensus/main' addParams( options: params.bcftools_consensus_options )
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PRE ; BCFTOOLS_INDEX as BCFTOOLS_INDEX_POST } from '../../modules/nf-core/modules/bcftools/index/main'     addParams( options: params.bcftools_index_options )
include { BCFTOOLS_VIEW }                                                                from '../../modules/nf-core/modules/bcftools/view/main'      addParams( options: params.bcftools_view_options )
include { FREEBAYES }                                                                    from '../../modules/nf-core/modules/freebayes/main'          addParams( options: params.freebayes_options )
include { PYDAMAGE_ANALYZE }                                                             from '../../modules/nf-core/modules/pydamage/analyze/main'   addParams( options: params.pydamage_analyze_options )
include { PYDAMAGE_FILTER }                                                              from '../../modules/nf-core/modules/pydamage/filter/main'    addParams( options: params.pydamage_filter_options )
include { SAMTOOLS_FAIDX as FAIDX}                                                       from '../../modules/nf-core/modules/samtools/faidx/main'     addParams( options: params.samtools_faidx_options )

workflow ANCIENT_DNA_ASSEMLY_VALIDATION {
    take:
        input //channel: [val(meta), path(contigs), path(bam), path(bam_index)]
    main:
        PYDAMAGE_ANALYZE(input.map {item -> [item[0], item[2], item[3]]})
        PYDAMAGE_FILTER(PYDAMAGE_ANALYZE.out.csv)
        FAIDX(input.map { item -> [ item[1] ] })
        FREEBAYES (input.map { item -> [item[0], item[2], item[3], [], []] }, input.map { item -> [ item[1] ] }, FAIDX.out.fai, [], [], [], [] )
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

