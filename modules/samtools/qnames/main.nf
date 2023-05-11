process QNAMES {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    } else {
        container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/samtools:1.14--hb421002_0' : null}"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"),  emit: qname
    path  "versions.yml",            emit: versions

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
        | uniq > ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}