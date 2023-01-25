process GENERATE_FILTERED_QNAMES {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.filtered_qnames.txt"), emit: qname

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def target_qnames = "${meta.id}_merged.target_qnames.txt"
    def filter_qnames = "${meta.id}_merged.filter_qnames.txt"
    """
    comm -2 -3 $target_qnames $filter_qnames > ${prefix}.filtered_qnames.txt
    """
}