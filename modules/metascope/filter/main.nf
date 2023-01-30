process METASCOPE_FILTER {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' ? 'file:///research/project/shared/benoukraf_lab/.singularity_cache/metascope.sif' : null}"

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