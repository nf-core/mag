//MODULES
include { MEGAHIT              } from '../../../modules/nf-core/megahit/main'
include { SPADES as METASPADES } from '../../../modules/nf-core/spades/main'

workflow SHORTREAD_ASSEMBLY {
    take:
    ch_short_reads_grouped // [val(meta), path(fastq1), path(fastq2)] (mandatory)
    ch_short_reads_spades // [val(meta), path(fastq1)]               (mandatory)

    main:
    ch_versions = channel.empty()
    ch_assembled_contigs = channel.empty()
    ch_publish = channel.empty()

    if (!params.single_end && !params.skip_spades) {
        METASPADES(ch_short_reads_spades.map { meta, reads -> [meta, reads, [], []] }, [], [])
        ch_spades_assemblies = (params.spades_downstreaminput == 'contigs' ? METASPADES.out.contigs : METASPADES.out.scaffolds).map { meta, assembly ->
            def meta_new = meta + [assembler: 'SPAdes']
            [meta_new, assembly]
        }
        ch_versions = ch_versions.mix(METASPADES.out.versions)

        ch_assembled_contigs = ch_assembled_contigs.mix(ch_spades_assemblies)
        ch_publish = ch_publish
            .mix(METASPADES.out.scaffolds, METASPADES.out.contigs, METASPADES.out.transcripts, METASPADES.out.gene_clusters, METASPADES.out.gfa, METASPADES.out.warnings, METASPADES.out.log)
            .map { meta, files -> [meta + [assembler: 'SPAdes'], files] }
    }

    if (!params.skip_megahit) {
        MEGAHIT(ch_short_reads_grouped)
        ch_megahit_assemblies = MEGAHIT.out.contigs.map { meta, assembly ->
            def meta_new = meta + [assembler: 'MEGAHIT']
            [meta_new, assembly]
        }
        ch_versions = ch_versions.mix(MEGAHIT.out.versions)

        ch_assembled_contigs = ch_assembled_contigs.mix(ch_megahit_assemblies)
        ch_publish = ch_publish.mix(
            MEGAHIT.out.contigs.mix(MEGAHIT.out.log).map { meta, files -> [meta + [assembler: 'MEGAHIT'], files] }
        )
    }

    emit:
    assembled_contigs = ch_assembled_contigs
    versions          = ch_versions
    publish           = ch_publish // [val(meta), path(files)]
}
