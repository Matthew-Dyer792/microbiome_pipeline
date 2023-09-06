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

include { METASCOPE_ID      } from '../modules/metascope/id/main'

process MERGE_METASCOPE_ID {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "anaconda::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/biopython:1.78' : null}"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'file:///research/project/shared/benoukraf_lab/.singularity_cache/biopython.sif' : '' }"

    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("${meta.id}.metascope_id.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def merged_metascope_id = "${meta.id}_merged.metascope_id.csv"
    def sorted_metascope_id = "${meta.id}.sorted.metascope_id.csv"
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *.csv > $merged_metascope_id

    (head -n 1 $merged_metascope_id && tail -n +2 $merged_metascope_id | sort -k1,1n -k3,3n) > $sorted_metascope_id

    python ${projectDir}/bin/merge.py \\
        -f $sorted_metascope_id \\
        -o ./
    """
}

process METASCOPE_TO_KRONA {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 1

    conda (params.enable_conda ? "anaconda::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/biopython:1.78' : null}"

    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/lineage.py \\
        -f $csv \\
        -o .
    """
}

process PLOT_KRONA {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "biconda::krona=2.81" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1' : null}"

    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.krona.html"), emit: krona

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    ktImportText \\
        -o ${prefix}.krona.html \\
        -n Microbiome \\
        $tsv
    """
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow IDENTIFY {

    take:
    bam // channel: val( meta ), path( bam )

    main:

    //
    // MODULE: Run Metascope's EM algorithm to identify the microbes
    //
    METASCOPE_ID (bam)

    // // refactor id to allow grouping of qname files by sample id
    // METASCOPE_ID.out.csv
    //     .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
    //     .groupTuple()
    //     .set{ bulk_ids }

    // //
    // // MODULE: Run Metascope's EM algorithm to identify the microbes
    // //
    // MERGE_METASCOPE_ID (bulk_ids)

    // //
    // // MODULE: Run Metascope's EM algorithm to identify the microbes
    // //
    // METASCOPE_TO_KRONA (MERGE_METASCOPE_ID.out.csv)

    // //
    // // MODULE: Run Metascope's EM algorithm to identify the microbes
    // //
    // PLOT_KRONA (METASCOPE_TO_KRONA.out.tsv)

    // nothing to emit as of now
}

/*
========================================================================================
    THE END
========================================================================================
*/