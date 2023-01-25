#!/usr/bin/env nextflow
/*
========================================================================================
    Microbiome Analysis Pipeline
========================================================================================
    take in cleaned up fastq of short-read or long-read format and produce 
    the quantity of every refseq stored bacteria and/or virus along with an
    optional Krona plot. Eventually we will add phyolgency contruction as well 
    hopefully leveraging CAIR GPU clusters
========================================================================================
    Github : TOO BE ADDED
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WorkflowMain.initialise(workflow, params, log)

// // Check if --input file is empty
// ch_input = file(params.input, checkIfExists: true)
// if (ch_input.isEmpty()) {exit 1, "File provided with --input is empty: ${ch_input.getName()}!"}

// Read in fastq from --input file
Channel
    .fromFilePairs(params.input, size: params.single_end ? 1 : 2, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .fastq files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }
    .map { it -> tuple([id: "${it[0]}"], it[1]) }
    .set { fastq_files }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// if (params.input_type == input_type) {
//     if (params.input_type == 'sra') {
//         include { SRA } from './workflows/sra'
//     } else if (params.input_type == 'synapse') {
//         include { SYNAPSE } from './workflows/synapse'
//     }
// } else {
//     exit 1, "Ids auto-detected as ${input_type}. Please provide '--input_type ${input_type}' as a parameter to the pipeline!"
// }

if (params.workflow == 'ont_long_reads') {
    include { ONT_LONG_READS } from './workflows/ont_long_reads'
// } else if (params.workflow == 'trimming') {
//     include { TRIMMING } from './workflows/trimming'
// } else if (params.workflow == 'align') {
//     include { TRIMMING } from './workflows/align'
}

include { IDENTIFY } from './workflows/identify'

//
// WORKFLOW: Run main methmotif pipeline depending on the step provided
//
workflow     MICROBIOME_PIPELINE {

    //
    // WORKFLOW: run the pre-alignment qc
    //
    if (params.workflow == 'ont_long_reads') {
        ONT_LONG_READS ( fastq_files )

        ONT_LONG_READS.out.id
            .set{ to_identify }

    // //
    // // WORKFLOW: trim the fastq files if necessary
    // //
    // } else if (params.workflow == 'trimming') {
    //     TRIMMING (  )

    // //
    // // WORKFLOW: trim the fastq files if necessary
    // //
    // } else if (params.workflow == 'align') {
    //     ALIGN (  )

    // to_identify.view()

    IDENTIFY ( to_identify )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    MICROBIOME_PIPELINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/