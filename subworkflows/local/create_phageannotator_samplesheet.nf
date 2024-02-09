include { CAT_CAT               } from '../../modules/nf-core/cat/cat/main'

workflow CREATE_PHAGEANNOTATOR_SAMPLESHEET {
    take:
        short_reads //channel: [val(meta), path(fastq_1), path(fastq_2)]
        assemblies  //channel: [val(meta), path(fasta)]
    main:
        ch_versions = Channel.empty()

        // combine assemblies by sample/group if multiple assembly methods were used
        ch_assemblies = assemblies
            .map {
                meta, fasta ->
                    def meta_new = meta.subMap('id')
                [ meta_new, fasta ]
            }
            .groupTuple()

        //
        // MODULE: Combine all assemblies from a sample into one FastA file
        //
        ch_combined_assemblies = CAT_CAT ( ch_assemblies ).file_out
        ch_versions = ch_versions.mix( CAT_CAT.out.versions )

        // if no coassembly, join FastQ and FastA by ID
        if ( !params.coassemble_group ){
            ch_combined_assemblies_remap = ch_combined_assemblies
                .map {
                    meta, fasta ->
                        def id          = meta.id

                        return [ id, fasta ]
                }

            short_reads
                .map {
                    meta, fastq ->
                        def id          = meta.id
                        def group       = meta.group
                        def single_end  = meta.single_end

                        return [ id, group, single_end, fastq ]
                }.join ( ch_combined_assemblies_remap )
                .map {
                    id, group, single_end, fastq, fasta ->
                        def reads   = fastq instanceof List ? fastq.flatten() : [ fastq ]
                        def meta    = [:]

                        meta.id         = id
                        meta.group      = group
                        meta.single_end = single_end
                        meta.fastq_1    = reads[0]
                        meta.fastq_2    = !meta.single_end ? reads[1] : ''
                        meta.fasta      = fasta ? fasta : ''

                        return meta
                }
                .set { ch_mag_metadata }
        } else {
            // if coassembly was used, join FastQ and FastA by group
            ch_combined_assemblies_remap = ch_combined_assemblies
                .map {
                    meta, fasta ->
                        def group = meta.id.split('group-')

                            return [ group[1], fasta ]
                }
            short_reads
                .map {
                    meta, fastq ->
                        def id          = meta.id
                        def group       = meta.group
                        def single_end  = meta.single_end

                            return [ group, id, single_end, fastq ]
                }
                .combine ( ch_combined_assemblies_remap, by:0 )
                .map {
                    group, id, single_end, fastq, fasta ->
                        def reads   = fastq instanceof List ? fastq.flatten() : [ fastq ]
                        def meta    = [:]

                        meta.id         = id
                        meta.group      = group
                        meta.single_end = single_end
                        meta.fastq_1    = reads[0]
                        meta.fastq_2    = !meta.single_end ? reads[1] : ''
                        meta.fasta      = fasta ? fasta : ''

                        return meta
                }
                .set { ch_mag_metadata }
        }

        // Create samplesheet for each sample using meta information
        ch_mag_id_samplesheets = ch_mag_metadata.collectFile() { meta ->
            [ "${meta.id}_phageannotator_samplesheet.csv", "sample,group,fastq_1,fastq_2,fasta" + '\n' + "${meta.id},${meta.group},${meta.fastq_1},${meta.fastq_2},${meta.fasta}" + '\n' ]
        }

        // Merge samplesheet across all samples for the pipeline
        ch_mag_id_samplesheets.collectFile(name: "phageannotator_samplesheet.csv", keepHeader:true, skip:1, storeDir:"${params.outdir}/downstream_samplesheets/")

    emit:
        versions            = ch_versions       // channel: [ versions.yml ]
}
