// SUBWORKFLOWS
include { SHORTREAD_ASSEMBLY } from './shortread_assembly'
include { LONGREAD_ASSEMBLY  } from './longread_assembly'
include { HYBRID_ASSEMBLY    } from './hybrid_assembly'

// MODULES
include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS          } from '../../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                     } from '../../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                  } from '../../modules/local/pool_single_reads'
include { GUNZIP as GUNZIP_SHORTREAD_ASSEMBLIES                 } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_LONGREAD_ASSEMBLIES                  } from '../../modules/nf-core/gunzip'

workflow ASSEMBLY {
    take:
    ch_short_reads           // [ [meta] , fastq1, fastq2] (mandatory)
    ch_long_reads            // [ [meta] , fastq] (mandatory)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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
                def meta = [:]
                meta.id = "group-${group}"
                meta.group = group
                meta.single_end = assemble_as_single
                if (assemble_as_single) {
                    [meta, reads.collect { it }, []]
                }
                else {
                    [meta, reads.collect { it[0] }, reads.collect { it[1] }]
                }
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
                [meta, reads.collect { it }]
            }
    }
    else {
        ch_short_reads_grouped = ch_short_reads
            .filter { it[0].single_end }
            .map { meta, reads -> [meta, [reads], []] }
            .mix(
                ch_short_reads.filter { !it[0].single_end }.map { meta, reads -> [meta, [reads[0]], [reads[1]]] }
            )
        ch_long_reads_grouped = ch_long_reads
    }

    if (!params.skip_spades || !params.skip_spadeshybrid) {
        if (params.coassemble_group) {
            if (params.bbnorm) {
                ch_short_reads_spades = ch_short_reads_grouped.map { [it[0], it[1]] }
            }
            else {
                POOL_SHORT_SINGLE_READS(
                    ch_short_reads_grouped.filter { it[0].single_end }
                )
                POOL_PAIRED_READS(
                    ch_short_reads_grouped.filter { !it[0].single_end }
                )
                ch_short_reads_spades = POOL_SHORT_SINGLE_READS.out.reads.mix(POOL_PAIRED_READS.out.reads)
            }
        }
        else {
            ch_short_reads_spades = ch_short_reads
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid) {
            POOL_LONG_READS(ch_long_reads_grouped)
            ch_long_reads_spades = POOL_LONG_READS.out.reads
        }
        else {
            ch_long_reads_spades = Channel.empty()
        }
    }
    else {
        ch_short_reads_spades = Channel.empty()
        ch_long_reads_spades = Channel.empty()
    }

    /*
    ================================================================================
                                        Assembly
    ================================================================================
    */

    ch_shortread_assembled_contigs = Channel.empty()
    ch_longread_assembled_contigs = Channel.empty()

    // SHORTREAD ASSEMBLY
    SHORTREAD_ASSEMBLY(
        ch_short_reads_grouped,
        ch_short_reads_spades
    )
    ch_versions = ch_versions.mix(SHORTREAD_ASSEMBLY.out.versions)

    // HYBRID ASSEMBLY
    HYBRID_ASSEMBLY(
        ch_short_reads_spades,
        ch_long_reads_spades
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
    longread_assemblies = ch_longread_assemblies
    versions = ch_versions
    multiqc_files = ch_multiqc_files

}
