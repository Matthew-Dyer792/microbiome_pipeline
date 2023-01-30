process METASCOPE_FILTER {
    tag "$meta.id"
    label 'process_high'
    
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file:///research/project/shared/benoukraf_lab/.singularity_cache/metascope.sif' : '' }"

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