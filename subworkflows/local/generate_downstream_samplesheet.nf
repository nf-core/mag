include { CAT_CAT               } from '../../modules/nf-core/cat/cat/main'

workflow GENERATE_DOWNSTREAM_SAMPLESHEET {
    take:
        downstream_nfcore_pipelines // val: [ nf-core-pipeline, OPTIONAL: other-nf-core-pipelines ]
        short_reads                 // channel: [val(meta), path(fastq_1), path(fastq_2)]
        assemblies                  // channel: [val(meta), path(fasta)]

    main:

        ch_versions = Channel.empty()

        //
        // Create a samplesheet for nf-core/phageannotator
        //
        if ( 'phageannotator' in downstream_nfcore_pipelines ) {

            if ( params.samplesheet_combine_assemblers ) {
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
            } else {
                ch_combined_assemblies = assemblies
            }

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
                    // Save reads and assemblies to outdir so that they are in a stable location
                    file(meta.fastq_1.toUriString(), checkIfExists: true).copyTo("${params.outdir}/downstream_samplesheets/fastq/${meta.fastq_1.name}")
                    file(meta.fasta, checkIfExists: true).copyTo("${params.outdir}/downstream_samplesheets/fasta/${meta.fasta.name}")
                    if ( !meta.single_end ){
                        file(meta.fastq_2.toUriString(), checkIfExists: true).copyTo("${params.outdir}/downstream_samplesheets/fastq/${meta.fastq_2.name}")
                        [ "${meta.id}_phageannotator_samplesheet.csv",
                            "sample,group,fastq_1,fastq_2,fasta" +
                            '\n' +
                            "${meta.id},${meta.group}," +
                            file("${params.outdir}/downstream_samplesheets/fastq/${meta.fastq_1.name}").toString() + "," +
                            file("${params.outdir}/downstream_samplesheets/fastq/${meta.fastq_2.name}").toString() + "," +
                            file("${params.outdir}/downstream_samplesheets/fasta/${meta.fasta.name}").toString() +
                            '\n'
                        ]
                    } else {
                        // Create samplesheet for each sample using meta information
                        [ "${meta.id}_phageannotator_samplesheet.csv",
                            "sample,group,fastq_1,fastq_2,fasta" +
                            '\n' +
                            "${meta.id},${meta.group}," +
                            file("${params.outdir}/downstream_samplesheets/fastq/${meta.fastq_1.name}").toString() + "," +
                            "," +
                            file("${params.outdir}/downstream_samplesheets/fasta/${meta.fasta.name}").toString() +
                            '\n'
                        ]
                    }
                }

            // Merge samplesheet across all samples for the pipeline
            ch_mag_id_samplesheets.collectFile(name: "phageannotator_samplesheet.csv", keepHeader:true, skip:1, storeDir:"${params.outdir}/downstream_samplesheets/")
        }

    emit:
        versions            = ch_versions       // channel: [ versions.yml ]
}
