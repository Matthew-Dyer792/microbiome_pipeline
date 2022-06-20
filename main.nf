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

// allow for either short reads (bowtie2) or long reads (minimap2)
if (!params.ont_long_reads) {
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
} else {
    // channel for the target genome index files
    Channel
        .fromPath(params.target_index, type: 'file', checkIfExists: true)
        .set{ target_index }

    // channel for composition of fastq reads and target indexed genomes
    fastq_files
        .combine(target_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            params.single_end ? file("${it[1][0]}") : [file("${it[1][0]}"), file("${it[1][1]}")],
            file("${it[2]}")) }
        .set{ align_files }

    // channel for the filter genome index files
    Channel
        .fromPath(params.filter_index, type: 'file', checkIfExists: true)
        .set { filter_index }
}




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
    tuple val(meta), path("*merged.target.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.target.bam"
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

process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    cpus 28
    memory '56 GB'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*bam"),              emit: bam
    tuple val(meta), path("*.minimap2.log"),    emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fastq   = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    def cpus    = task.cpus - 1
    """
    minimap2 \\
        $args \\
        -t $cpus \\
        -a $index \\
        $fastq \\
        2> ${prefix}.minimap2.log \\
        | samtools view -@ $cpus $args2 -bhS -o ${prefix}.bam
    """
}

process BAM_TO_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    publishDir "${params.outdir}/${meta.id}_results", mode: 'copy'

    cpus 16
    memory '32 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2    = task.ext.args2 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.target.fastq"
    """
    samtools fastq \\
        $args \\
        -@ $task.cpus \\
        $bam \\
        > $output
    
    pigz \\
        $args2 \\
        -p $task.cpus\\
        *.fastq
    """
}

process TARGET_QNAMES {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.target_qnames.txt"), emit: qname

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    samtools view \\
        -@ $task.cpus \\
        $args \\
        $bam \\
        | cut -f1 \\
        | sort -T . \\
        | uniq > ${prefix}.target_qnames.txt
    """
}

process MINIMAP2_FILTER {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    cpus 28
    memory '56 GB'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.txt"),             emit: qname
    tuple val(meta), path("*.minimap2.log"),    emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fastq   = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    def cpus    = task.cpus - 1
    """
    minimap2 \\
        $args \\
        -t $cpus \\
        -a $index \\
        $fastq \\
        2> ${prefix}.minimap2.log \\
        | samtools view -@ $cpus $args2 \\
        | cut -f1 \\
        | sort -T . \\
        | uniq > ${prefix}.txt
    """
}

process QNAMES_MERGE {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    cpus 1
    memory '8 GB'

    input:
    tuple val(meta), path(qnames)

    output:
    tuple val(meta), path("*.filter_qnames.txt"), emit: qname

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    cat $qnames \\
        | sort -T . \\
        | uniq > ${prefix}.filter_qnames.txt
    """
}

process FILTER_BAM {
    tag "$meta.id"
    label 'process_medium'

    // errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    // maxRetries 3

    errorStrategy 'finish'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    publishDir "${params.outdir}/${meta.id}_results", mode: 'copy'

    cpus 16
    memory '32 GB'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def target_qnames = "${meta.id}.target_qnames.txt"
    def filter_qnames = "${meta.id}_merged.filter_qnames.txt"
    def filtered_qnames = "${meta.id}.filtered_qnames.txt"
    def bam = "${meta.id}_merged.target.bam"
    """
    comm -2 -3 $target_qnames $filter_qnames > $filtered_qnames

    samtools view -@ $task.cpus -bh -N $filtered_qnames -o ${prefix}.filtered.bam $bam
    """
}

process METASCOPE_TO_KRONA {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/pathoscope/metascope/minimap2-2.24.sif' : '' }"

    publishDir "${params.outdir}/${meta.id}_results", mode: 'copy'

    cpus 1
    memory '8 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2    = task.ext.args2 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.target.fastq"
    """
    samtools fastq \\
        $args \\
        -@ $task.cpus \\
        $bam \\
        > $output
    
    pigz \\
        $args2 \\
        -p $task.cpus\\
        *.fastq
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
    if (!params.ont_long_reads) {
        METASCOPE_ALIGN (align_files)

        // refactor id to allow grouping of bam files by sample id
        METASCOPE_ALIGN.out.bam
            .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
            .groupTuple()
            .set{ aligned_bam }
    } else {
        MINIMAP2_ALIGN (align_files)

        // refactor id to allow grouping of bam files by sample id
        MINIMAP2_ALIGN.out.bam
            .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
            .groupTuple()
            .set{ aligned_bam }
    }
    
    BAM_MERGE (aligned_bam)

    if (!params.ont_long_reads) {
        // add cartesian mulitply the sample fastqs to the indexed target genomes
        METASCOPE_ALIGN.out.bam
            .combine(filter_index)
            .map{ it-> 
                tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
                file("${it[1]}"),
                [file("${it[2]}"), file("${it[3]}"), file("${it[4]}"), file("${it[5]}"), file("${it[6]}"), file("${it[7]}")]) }
            .set{ filter_files }

        filter_files.view()

        METASCOPE_FILTER (filter_files)

        METASCOPE_FILTER.out.bam
            .set{ id }
    } else {
        TARGET_QNAMES (BAM_MERGE.out.bam)
        BAM_TO_FASTQ (BAM_MERGE.out.bam)

        // cartesian mulitply the sample fastqs to the indexed filter genomes
        BAM_TO_FASTQ.out.fastq
            .combine(filter_index)
            .map{ it-> 
                tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
                file("${it[1]}"),
                file("${it[2]}")) }
            .set{ filter_files }

        MINIMAP2_FILTER (filter_files)

        // refactor id to allow grouping of qname files by sample id
        MINIMAP2_FILTER.out.qname
            .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
            .groupTuple()
            .set{ filter_qname }

        QNAMES_MERGE (filter_qname)

        BAM_MERGE.out.bam
            .mix( TARGET_QNAMES.out.qname, QNAMES_MERGE.out.qname )
            .groupTuple()
            .set{ to_filter }

        FILTER_BAM (to_filter)

        FILTER_BAM.out.bam
            .set{ id }
    }

    METASCOPE_ID (id)
    // METASCOPE_TO_KRONA ()
    // PLOT_KRONA
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
