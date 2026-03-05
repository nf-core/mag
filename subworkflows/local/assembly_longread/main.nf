/*
* LONGREAD_ASSEMBLY: Assembly of long reads
*/

include { FLYE         } from '../../../modules/nf-core/flye/main'
include { METAMDBG_ASM } from '../../../modules/nf-core/metamdbg/asm/main'

workflow LONGREAD_ASSEMBLY {
    take:
    ch_long_reads // [val(meta), path(fastq)] (mandatory)

    main:
    ch_assembled_contigs = channel.empty()
    ch_versions = channel.empty()

    if (!params.skip_flye) {

        ch_long_reads_flye_input = ch_long_reads.multiMap { meta, reads ->
            def fly_mode = ""
            if (meta.lr_platform == "OXFORD_NANOPORE") {
                fly_mode = "--nano-raw"
            }
            else if (meta.lr_platform == "OXFORD_NANOPORE_HQ") {
                fly_mode = "--nano-hq"
            }
            else if (meta.lr_platform == "PACBIO_HIFI") {
                fly_mode = "--pacbio-hifi"
            }
            else if (meta.lr_platform == "PACBIO_CLR") {
                fly_mode = "--pacbio-raw"
            }
            else {
                log.error("[nf-core/mag]: ERROR - unknown lr_platform provided to Flye!")
            }
            reads: [meta, reads]
            mode: fly_mode
        }

        FLYE(
            ch_long_reads_flye_input.reads,
            ch_long_reads_flye_input.mode,
        )
        ch_versions = ch_versions.mix(FLYE.out.versions)

        ch_flye_assemblies = FLYE.out.fasta.map { meta, assembly ->
            def meta_new = meta + [assembler: "FLYE"]
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_flye_assemblies)
    }

    if (!params.skip_metamdbg) {

        ch_long_reads_metamdbg_input = ch_long_reads.multiMap { meta, reads ->
            def metamdbg_mode = ""
            if (meta.lr_platform == "OXFORD_NANOPORE") {
                metamdbg_mode = "ont"
            }
            else if (meta.lr_platform == "OXFORD_NANOPORE_HQ") {
                metamdbg_mode = "ont"
            }
            else if (meta.lr_platform == "PACBIO_HIFI") {
                metamdbg_mode = "hifi"
            }
            else if (meta.lr_platform == "PACBIO_CLR") {
                metamdbg_mode = "hifi"
            }
            else {
                log.error("[nf-core/mag]: ERROR - unknown lr_platform provided to MetaMDBG!")
            }
            reads: [meta, reads]
            mode: metamdbg_mode
        }

        METAMDBG_ASM(
            ch_long_reads_metamdbg_input.reads,
            ch_long_reads_metamdbg_input.mode,
        )
        ch_versions = ch_versions.mix(METAMDBG_ASM.out.versions)

        ch_metamdbg_assemblies = METAMDBG_ASM.out.contigs.map { meta, assembly ->
            def meta_new = meta + [assembler: "METAMDBG"]
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_metamdbg_assemblies)
    }

    emit:
    assembled_contigs = ch_assembled_contigs
    versions          = ch_versions
}
