process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        conda "bioconda::samtools=1.14"
    } else {
        container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/samtools:1.14--hb421002_0' : null}"
    }

    input:
    tuple val(meta), path(bam)
    val(interleave)

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz"),          emit: fastq,        optional:true
    tuple val(meta), path("*_interleaved.fastq.gz"),    emit: interleaved,  optional:true
    tuple val(meta), path("*_singleton.fastq.gz"),      emit: singleton,    optional:true
    tuple val(meta), path("*_other.fastq.gz"),          emit: other,        optional:true
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def output = ( interleave && ! meta.single_end ) ? "> ${prefix}_interleaved.fastq.gz" :
        meta.single_end ? "-1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz" :
        "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
    """
    samtools fastq \\
        $args \\
        -@ $task.cpus \\
        -0 ${prefix}_other.fastq.gz \\
        $bam \\
        $output
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}