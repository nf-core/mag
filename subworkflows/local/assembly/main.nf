// SUBWORKFLOWS
include { SHORTREAD_ASSEMBLY                    } from '../assembly_shortread/main'
include { LONGREAD_ASSEMBLY                     } from '../assembly_longread/main'
include { HYBRID_ASSEMBLY                       } from '../assembly_hybrid/main'

// MODULES
include { CAT_FASTQ as POOL_SHORT_READS         } from '../../../modules/nf-core/cat/fastq'
include { CAT_FASTQ as POOL_LONG_READS          } from '../../../modules/nf-core/cat/fastq'
include { GUNZIP as GUNZIP_SHORTREAD_ASSEMBLIES } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_LONGREAD_ASSEMBLIES  } from '../../../modules/nf-core/gunzip'

workflow ASSEMBLY {
    take:
    ch_short_reads // [val(meta), path(fastq1), path(fastq2)] (mandatory)
    ch_long_reads  // [val(meta), path(fastq)]                (mandatory)

    main:

    ch_versions = channel.empty()

    /*
    ================================================================================
                                    Assembly preparation
    ================================================================================
    */

    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads_grouped = ch_short_reads
            .map { meta, reads -> [meta.group, meta, reads] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def assemble_as_single = params.single_end || (params.bbnorm && params.coassemble_group)
                def meta = [
                    id: "group-${group}",
                    group: group,
                    single_end: assemble_as_single,
                    sr_platform: metas.sr_platform[0]
                ]
                if (assemble_as_single) {
                    [meta, reads.sort { files -> files[0].getName() }, []]
                }
                else {
                    [meta] + reads.sort { files -> files[0].getName() }.transpose()
                }
            }

        // We have to merge reads together to match tuple structure of POOL_SHORT_READS/
        // This MUST be in a interleaved structure (s1_r1, s1_r2, s2_r1, s2_r2, ...)
        // So we zip the lists of R1s and R2s, consistent order is ensured by the existing structure
        ch_short_reads_grouped_for_pooling = ch_short_reads_grouped.map { meta, reads1, reads2 ->
            [meta, [reads1, reads2].transpose().flatten()]
        }

        // long reads
        // group and set group as new id
        ch_long_reads_grouped = ch_long_reads
            .map { meta, reads -> [meta.group, meta, reads] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def meta = [:]
                meta.id = "group-${group}"
                meta.group = group
                meta.lr_platform = metas.lr_platform[0]
                [meta, reads]
            }
    }
    else {
        ch_short_reads_grouped = ch_short_reads.map { meta, reads ->
            if (meta.single_end) {
                [meta, [reads], []]
            }
            else {
                [meta, [reads[0]], [reads[1]]]
            }
        }

        ch_long_reads_grouped = ch_long_reads
    }

    if (!params.skip_spades || !params.skip_spadeshybrid) {
        if (params.coassemble_group) {
            if (params.bbnorm) {
                // When doing co-assembly and using bbnorm, all sample reads get pooled in a single file
                // That's why we can drop R2 here (it's empty)
                ch_short_reads_spades = ch_short_reads_grouped.map { meta, r1, _r2 -> [meta, r1] }
            }
            else {
                POOL_SHORT_READS(ch_short_reads_grouped_for_pooling)
                ch_versions = ch_versions.mix(POOL_SHORT_READS.out.versions)
                ch_short_reads_spades = POOL_SHORT_READS.out.reads
            }
        }
        else {
            ch_short_reads_spades = ch_short_reads
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid) {

            ch_long_reads_grouped_for_pool = ch_long_reads_grouped
                .map { meta, reads -> [meta.id, meta, reads] }
                .combine(ch_short_reads_grouped.map { meta, _reads1, _reads2 -> [meta.id] }, by: 0)
                .map { _id, lr_meta, lr_reads -> [lr_meta, lr_reads] }
            //make sure no long reads are pooled for spades if there are no short reads

            POOL_LONG_READS(ch_long_reads_grouped_for_pool)
            ch_versions = ch_versions.mix(POOL_LONG_READS.out.versions)
            ch_long_reads_spades = POOL_LONG_READS.out.reads
        }
        else {
            ch_long_reads_spades = channel.empty()
        }
    }
    else {
        ch_short_reads_spades = channel.empty()
        ch_long_reads_spades = channel.empty()
    }

    /*
    ================================================================================
                                        Assembly
    ================================================================================
    */

    ch_shortread_assembled_contigs = channel.empty()
    ch_longread_assembled_contigs = channel.empty()

    // SHORTREAD ASSEMBLY
    SHORTREAD_ASSEMBLY(
        ch_short_reads_grouped,
        ch_short_reads_spades,
    )
    ch_versions = ch_versions.mix(SHORTREAD_ASSEMBLY.out.versions)

    // HYBRID ASSEMBLY
    HYBRID_ASSEMBLY(
        ch_short_reads_spades,
        ch_long_reads_spades,
    )
    ch_versions = ch_versions.mix(HYBRID_ASSEMBLY.out.versions)

    // LONGREAD ASSEMBLY
    LONGREAD_ASSEMBLY(
        ch_long_reads_grouped
    )
    ch_versions = ch_versions.mix(LONGREAD_ASSEMBLY.out.versions)

    ch_shortread_assembled_contigs = SHORTREAD_ASSEMBLY.out.assembled_contigs.mix(HYBRID_ASSEMBLY.out.assembled_contigs)
    ch_longread_assembled_contigs = LONGREAD_ASSEMBLY.out.assembled_contigs

    GUNZIP_SHORTREAD_ASSEMBLIES(ch_shortread_assembled_contigs)
    ch_versions = ch_versions.mix(GUNZIP_SHORTREAD_ASSEMBLIES.out.versions)
    ch_shortread_assemblies = GUNZIP_SHORTREAD_ASSEMBLIES.out.gunzip

    GUNZIP_LONGREAD_ASSEMBLIES(ch_longread_assembled_contigs)
    ch_versions = ch_versions.mix(GUNZIP_LONGREAD_ASSEMBLIES.out.versions)
    ch_longread_assemblies = GUNZIP_LONGREAD_ASSEMBLIES.out.gunzip

    emit:
    shortread_assemblies = ch_shortread_assemblies
    longread_assemblies  = ch_longread_assemblies
    versions             = ch_versions
}
