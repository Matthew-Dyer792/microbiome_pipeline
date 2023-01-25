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

include { SAMTOOLS_MERGE as FILTERED_BAM_MERGE  } from '../modules/samtools/merge/main'
include { MERGE_QNAMES as MERGE_TARGET_QNAMES   } from '../modules/command_line/merge_qnames/main'
include { MERGE_QNAMES as MERGE_FILTER_QNAMES   } from '../modules/command_line/merge_qnames/main'
include { GENERATE_FILTERED_QNAMES              } from '../modules/command_line/compare_qnames/main'
// include { BISMARK_ALIGN      } from '../modules/bismark/align'
// maybe add the samtools subworkflow

process MINIMAP2_ALIGN {
    tag "$meta.id"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "bioconda::minimap2=2.24 bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' : null}"

    // publishDir "${params.outdir}/${meta.id}_results/bam", mode: 'copy'

    cpus 12
    memory '24 GB'

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
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        -a $index \\
        $fastq \\
        2> ${prefix}.minimap2.log \\
        | samtools view -@ $task.cpus $args2 -bh -o ${prefix}.bam
    """
}

process MINIMAP2_FILTER {
    tag "$meta.id"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "bioconda::minimap2=2.24 bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' : null}"

    cpus 12
    memory '24 GB'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.txt"),             emit: qname
    tuple val(meta), path("*.minimap2.log"),    emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fastq   = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        -a $index \\
        $fastq \\
        2> ${prefix}.minimap2.log \\
        | samtools view -@ $task.cpus \\
        | cut -f1 \\
        | sort -T . \\
        | uniq > ${prefix}.txt
    """
}

process TARGET_QNAMES {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/samtools:1.14--hb421002_0' : null}"

    cpus 4
    memory '8 GB'

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

process FILTER_BAM {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/samtools:1.14--hb421002_0' : null}"

    // publishDir "${params.outdir}/${meta.id}_results", mode: 'copy'

    cpus 8
    memory '24 GB'

    input:
    tuple val(meta), path(bam), path(filtered_qnames)

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    mkdir test/

    samtools view -@ $task.cpus -bh -N $filtered_qnames -o ${prefix}.filtered.bam $bam
    """
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow ONT_LONG_READS {

    take:
    fastq_files // channel: val( meta ), path([ fastq ])

    main:

    //
    // Setup the index and alignment channels // move this to separate validation file
    //

    // channel for the target genome index files
    Channel
        .fromPath(params.target_index, type: 'file', checkIfExists: true)
        .set{ target_index }

    // channel for composition of fastq reads and target indexed genomes
    fastq_files
        .combine(target_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", target_index_prefix: "${it[0].index_prefix}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            params.single_end ? file("${it[1][0]}") : [file("${it[1][0]}"), file("${it[1][1]}")],
            file("${it[2]}")) }
        .set{ target_align_files }

    // target_align_files
    //     .first()
    //     .view()

    // channel for the filter genome index files
    Channel
        .fromPath(params.filter_index, type: 'file', checkIfExists: true)
        .set { filter_index }

    // channel for composition of fastq reads and filter indexed genomes
    fastq_files
        .combine(filter_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", target_index_prefix: "${it[0].index_prefix}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            params.single_end ? file("${it[1][0]}") : [file("${it[1][0]}"), file("${it[1][1]}")],
            file("${it[2]}")) }
        .set{ filter_align_files }

    // filter_align_files
    //     .first()
    //     .view()

    //
    // MODULE: Align the raw fastq files against the given microbe index files
    //
    MINIMAP2_ALIGN (target_align_files)

    //
    // MODULE: Extract the qnames from the target aligned reads
    //
    TARGET_QNAMES (MINIMAP2_ALIGN.out.bam)

    // refactor id to allow grouping of qname files by sample id
    TARGET_QNAMES.out.qname
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ target_qname }

    // target_qname.view()

    //
    // MODULE: Merge the target qnames and then return a file containing only the unique target qnames
    //
    MERGE_TARGET_QNAMES (target_qname)

    //
    // MODULE: Align the raw fastq files against the given filter index files (human, etc...)
    //
    MINIMAP2_FILTER (filter_align_files)

    // refactor id to allow grouping of bam files by sample id
    MINIMAP2_FILTER.out.qname
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ filter_qnames }

    // filter_qnames.view() // this needs to be merged as there maybe more than one filter genome

    //
    // MODULE: Merge the filter qnames and then return a file containing only the unique filter qnames
    //
    MERGE_FILTER_QNAMES (filter_qnames)

    // refactor id to allow grouping of qname files by sample id
    MERGE_TARGET_QNAMES.out.qname
        .mix( MERGE_FILTER_QNAMES.out.qname )
        .groupTuple()
        .set{ filtered_qname }

    // filtered_qname.view()

    //
    // MODULE: Merge the target & filter qnames together and then return a file containing only the qnames present in both
    //
    GENERATE_FILTERED_QNAMES (filtered_qname)

    // need to remap to set the 2 index of the tuple to a shared ID
    MINIMAP2_ALIGN.out.bam
        .map{ it->
            tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[0].index_prefix}"],
            file("${it[1]}"),
            "${it[0].id}") }
        .set{ target_bam }

    // need to remap to set the 2 index of the tuple to a shared ID
    GENERATE_FILTERED_QNAMES.out.qname
        .map{ it->
            tuple([id: "${it[0].id}"],
            file("${it[1]}"),
            "${it[0].id}") }
        .set{ filtered_qnames }

    // cartesian mulitply the aligned bams to the filtered qnames
    target_bam
        .combine( filtered_qnames, by: 2 ) // 2 is the index of the shared ID
        .map{ it-> 
            tuple([id: "${it[1].id}", single_end: params.single_end, index_prefix: "${it[1].index_prefix}"],
            file("${it[2]}"),
            file("${it[4]}")) }
        .set{ to_filter }

    // to_filter.first().view()

    //
    // MODULE: Filter the target bam file to remove all reads that aligned to the filter genomes
    //
    FILTER_BAM (to_filter)

    // refactor id to allow grouping of filtered bam files by sample id
    FILTER_BAM.out.bam
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ to_merge }

    // to_merge.view()

    //
    // MODULE: Merge the filtered bam files together
    //
    FILTERED_BAM_MERGE (to_merge)

    emit:
    id          = FILTERED_BAM_MERGE.out.bam      // channel: [ val(meta), bam   ]

}

/*
========================================================================================
    THE END
========================================================================================
*/