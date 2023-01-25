process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/samtools:1.14--hb421002_0' : null}"

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