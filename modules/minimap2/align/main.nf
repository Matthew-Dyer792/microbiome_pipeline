process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    if (params.enable_conda) {
        conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"
    } else {
        container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' : null}"
    }

    input:
    tuple val(meta), path(reads), path(index)
    val bam_format

    output:
    tuple val(meta), path("*.bam"),             optional: true, emit: bam
    tuple val(meta), path("*.paf"),             optional: true, emit: paf
    tuple val(meta), path("*.minimap2.log"),    emit: log
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fastq   = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    def bam_output = bam_format ? "-a | samtools sort | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $index \\
        $fastq \\
        2> ${prefix}.minimap2.log \\
        $bam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}