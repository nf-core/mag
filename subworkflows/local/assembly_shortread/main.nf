//MODULES
include { MEGAHIT              } from '../../../modules/nf-core/megahit/main'
include { SPADES as METASPADES } from '../../../modules/nf-core/spades/main'

workflow SHORTREAD_ASSEMBLY {
    take:
    ch_short_reads_grouped // [ [meta] , fastq1, fastq2] (mandatory)
    ch_short_reads_spades

    main:
    ch_versions = Channel.empty()
    ch_assembled_contigs = Channel.empty()

    if (!params.single_end && !params.skip_spades) {
        METASPADES(ch_short_reads_spades.map { meta, reads -> [meta, reads, [], []] }, [], [])
        ch_spades_assemblies = METASPADES.out.scaffolds.map { meta, assembly ->
            def meta_new = meta + [assembler: 'SPAdes']
            [meta_new, assembly]
        }
        ch_versions = ch_versions.mix(METASPADES.out.versions)

        ch_assembled_contigs = ch_assembled_contigs.mix(ch_spades_assemblies)
    }

    if (!params.skip_megahit) {
        MEGAHIT(ch_short_reads_grouped)
        ch_megahit_assemblies = MEGAHIT.out.contigs.map { meta, assembly ->
            def meta_new = meta + [assembler: 'MEGAHIT']
            [meta_new, assembly]
        }
        ch_versions = ch_versions.mix(MEGAHIT.out.versions)

        ch_assembled_contigs = ch_assembled_contigs.mix(ch_megahit_assemblies)
    }

    emit:
    assembled_contigs = ch_assembled_contigs
    versions          = ch_versions
}
