/*
* LONGREAD_ASSEMBLY: Assembly of long reads
*/

include { FLYE         } from '../../../modules/nf-core/flye/main'
include { METAMDBG_ASM } from '../../../modules/nf-core/metamdbg/asm/main'

workflow LONGREAD_ASSEMBLY {
    take:
    ch_long_reads // [ [meta] , fastq] (mandatory)

    main:
    ch_assembled_contigs = Channel.empty()
    ch_versions = Channel.empty()

    if (!params.skip_flye) {

        ch_long_reads_flye_input = ch_long_reads.multiMap { meta, reads ->
            reads: [meta, reads]
            mode: meta.lr_platform == "OXFORD_NANOPORE"
                ? "--nano-raw"
                : meta.lr_platform == "OXFORD_NANOPORE_HQ"
                    ? "--nano-hq"
                    : meta.lr_platform == "PACBIO_HIFI"
                        ? "--pacbio-hifi"
                        : meta.lr_platform == "PACBIO_CLR" ? "--pacbio-raw" : []
        }

        FLYE(
            ch_long_reads_flye_input.reads,
            ch_long_reads_flye_input.mode,
        )
        ch_versions = ch_versions.mix(FLYE.out)

        ch_flye_assemblies = FLYE.out.fasta.map { meta, assembly ->
            def meta_new = meta + [assembler: "FLYE"]
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_flye_assemblies)
    }

    if (!params.skip_metamdbg) {

        ch_long_reads_metamdbg_input = ch_long_reads.multiMap { meta, reads ->
            reads: [meta, reads]
            mode: meta.lr_platform == "OXFORD_NANOPORE"
                ? "ont"
                : meta.lr_platform == "OXFORD_NANOPORE_HQ"
                    ? "ont"
                    : meta.lr_platform == "PACBIO_HIFI"
                        ? "hifi"
                        : meta.lr_platform == "PACBIO_CLR" ? "hifi" : []
        }

        METAMDBG_ASM(
            ch_long_reads_metamdbg_input.reads,
            ch_long_reads_metamdbg_input.mode,
        )
        ch_versions = ch_versions.mix(METAMDBG_ASM.out)

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
