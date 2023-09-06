process METASCOPE_ID {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-metascope=1.0"
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
    // Default key from ACC account
    def ncbi_key  = params.ncbi_key ? params.ncbi_key : "522cc5cd1ca59343ba9f55282d9ecfe6c009"
    """
    Rscript --vanilla ${projectDir}/bin/MetaScope_id.R \\
        $ncbi_key \\
        $bam \\
        $aligner

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaScope: \$( cat sessionInfo.txt )
    END_VERSIONS
    """
}