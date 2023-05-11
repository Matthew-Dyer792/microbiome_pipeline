process SPLIT_BAM {
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' : null}"

    input:
    tuple val(meta), path(bam)
    val mem_limit from params.mem_limit // review this

    output:
    tuple val(meta), path("split_bams/*.bam"), emit: split_bams

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def bam_size = bam.size().div(1024) // get the size of the BAM file in kilobytes
    def n_splits = Math.ceil(bam_size.div(mem_limit)) // calculate the number of splits based on memory limit
    """
    # create output directory
    mkdir -p split_bams

    # split the bam file using Picard
    java -jar /usr/local/share/picard-3.0.0-1/picard.jar SplitSamByNumberOfReads \
        I=$bam \
        OUTPUT=split_bams \
        --OUT_PREFIX=$prefix \
        SPLIT_TO_N_FILES=$n_splits
    """
}