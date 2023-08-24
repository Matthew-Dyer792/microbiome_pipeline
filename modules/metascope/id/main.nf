process METASCOPE_ID {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' ? 'matthewdyer/metascope:1.0.0' : null}"

    input:
    tuple val(meta), path(bam, stageAs: 'bam/*')

    output:
    tuple val(meta), path("bam/*.csv"),     emit: csv
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def aligner   = params.ont_long_reads ? "other" : "bowtie2"
    """
    Rscript --vanilla ${projectDir}/bin/MetaScope_id.R \\
        $bam \\
        $aligner

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaScope: \$( cat sessionInfo.txt )
    END_VERSIONS
    """
}