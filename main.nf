#!/usr/bin/env nextflow
/*
========================================================================================
    methylKit output to CpA methylation.bed
========================================================================================
    take in methylKit CHH context text files and produce the trinucleotide context while
    keeping only CpA's using pysam and awk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE
========================================================================================
*/



/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// Create a fastq channel for input files
Channel
    .fromFilePairs(params.input, size: params.single_end ? 1 : 2, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .fastq files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }
    .map { it -> tuple([id: "${it[0]}"], it[1]) }
    .set { fastq_files }

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

/*
========================================================================================
    PIPELINE STEPS
========================================================================================
*/

process METASCOPE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/metascope.sif' : '' }"

    cpus 28
    memory '56 GB'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fastq   = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    def index_prefix  = "${meta.index_prefix}"
    """
    Rscript --vanilla ${projectDir}/bin/MetaScope_align.R \\
        $fastq \\
        $index_prefix \\
        $prefix \\
        $task.cpus
    """
}

process BAM_MERGE {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/metascope.sif' : '' }"

    publishDir "${params.outdir}/${meta.id}_results", mode: 'copy'

    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*merged.${meta.step}.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.bam"
    """
    samtools merge \\
        $args \\
        -@ $task.cpus \\
        -o $output \\
        $bam
    """
}

process METASCOPE_FILTER {
    tag "$meta.id"
    label 'process_high'
    
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/metascope.sif' : '' }"

    cpus 28
    memory '56 GB'

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_prefix  = "${meta.index_prefix}"
    """
    Rscript --vanilla ${projectDir}/bin/MetaScope_filter.R \\
        $bam \\
        $index_prefix \\
        $task.cpus
    """
}

process METASCOPE_ID {
    tag "$meta.id"
    label 'process_low'

    errorStrategy 'finish'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/metascope.sif' : '' }"

    publishDir "${params.outdir}/${meta.id}_results", mode: 'copy'

    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MetaScope_id.R \\
        $bam \\
        bowtie
    """
}

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

//
// WORKFLOW: Run main metascope analysis pipeline
//
workflow metascope_pipeline {
    METASCOPE_ALIGN (align_files)

    // refactor id to allow grouping of bam files by sample id
    METASCOPE_ALIGN.out.bam
        .map{ it -> tuple([id: "${it[0].id}", step: "target"], it[1]) }
        .groupTuple()
        .set{ aligned_bam }

    BAM_MERGE (aligned_bam)

    // add cartesian mulitply the sample fastqs to the indexed target genomes
    METASCOPE_ALIGN.out.bam
        .combine(filter_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            file("${it[1]}"),
            [file("${it[2]}"), file("${it[3]}"), file("${it[4]}"), file("${it[5]}"), file("${it[6]}"), file("${it[7]}")]) }
        .set{ filter_files }

    METASCOPE_FILTER (filter_files)

    // refactor id to allow grouping of bam files by sample id
    METASCOPE_FILTER.out.bam
        .map{ it -> tuple([id: "${it[0].id}", step: "filtered"], it[1]) }
        .groupTuple()
        .set{ filtered_bam }

    BAM_MERGE (filtered_bam)
    METASCOPE_ID (BAM_MERGE.out.bam)
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    metascope_pipeline ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
