process METASCOPE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' ? 'file:///research/project/shared/benoukraf_lab/.singularity_cache/metascope.sif' : null}"

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