// using output folder to keep results from inputs
process MERGE_QNAMES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(qnames)

    output:
    tuple val(meta), path("output/${task.ext.prefix}.txt"), emit: qname

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    mkdir output
    
    cat $qnames \\
        | sort -T . \\
        | uniq > output/${prefix}.txt
    """
}
