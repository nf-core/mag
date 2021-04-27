#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mag
========================================================================================

nf-core/mag Analysis Pipeline. Started 2018-05-22.
#### Homepage / Documentation
https://github.com/nf-core/mag
#### Authors
Hadrien Gourlé HadrienG <hadrien.gourle@slu.se> - hadriengourle.com>
Daniel Straub <d4straub@gmail.com>
Sabrina Krakau <sabrinakrakau@gmail.com>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/mag --input 'samplesheet.csv' -profile docker"
    // nextflow run nf-core/mag --input '*_R{1,2}.fastq.gz' -profile docker
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////

workflow {
    include { MAG } from './workflows/mag' addParams( summary_params: summary_params )
    MAG ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////