process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

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