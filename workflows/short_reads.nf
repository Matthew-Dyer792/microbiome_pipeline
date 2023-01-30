/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def valid_params = [
//     ena_metadata_fields : ['run_accession', 'experiment_accession', 'library_layout', 'fastq_ftp', 'fastq_md5']
// ]

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowSra.initialise(params, log, valid_params)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { SAMTOOLS_MERGE as BAM_MERGE           } from '../modules/samtools/merge/main'
// include { BISMARK_ALIGN      } from '../modules/bismark/align'
// maybe add the samtools subworkflow

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SHORT_READS {

    take:
    fastq_files // channel: val( meta ), path([ fastq ])

    main:

    //
    // Setup the index and alignment channels // move this to separate validation file
    //

    // channel for the target genome index files
    Channel
        .fromPath(params.target_index, type: 'file', checkIfExists: true)
        .toSortedList()
        .flatten()
        .collate( 6 )
        .set{ target_index }

    // channel for composition of fastq reads and target indexed genomes
    fastq_files
        .combine(target_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            params.single_end ? file("${it[1][0]}") : [file("${it[1][0]}"), file("${it[1][1]}")],
            [file("${it[2]}"), file("${it[3]}"), file("${it[4]}"), file("${it[5]}"), file("${it[6]}"), file("${it[7]}")]) }
        .set{ align_files }

    // channel for the filter genome index files
    Channel
        .fromPath(params.filter_index, type: 'file', checkIfExists: true)
        .toSortedList()
        .flatten()
        .collate( 6 )
        .set { filter_index }

    //
    // MODULE: Align the raw fastq files against the given microbe index files
    //
    METASCOPE_ALIGN (align_files)

    // refactor id to allow grouping of bam files by sample id
    METASCOPE_ALIGN.out.bam
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ aligned_bam }

    //
    // MODULE: Merge the filtered bam files together
    //
    BAM_MERGE (aligned_bam)


    // this needs a rework to multiply by id 
    // add cartesian mulitply the sample fastqs to the indexed target genomes
    METASCOPE_ALIGN.out.bam
        .combine(filter_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            file("${it[1]}"),
            [file("${it[2]}"), file("${it[3]}"), file("${it[4]}"), file("${it[5]}"), file("${it[6]}"), file("${it[7]}")]) }
        .set{ filter_files }

    filter_files.view()

    //
    // MODULE: Align the new fastq files against the given filter index files (human, etc...) and remove all reads that aligned to both genomes
    //
    METASCOPE_FILTER (filter_files)

    emit:
    id          = METASCOPE_FILTER.out.bam      // channel: [ val(meta), bam   ]

}

/*
========================================================================================
    THE END
========================================================================================
*/