process METASCOPE_ID {
    tag "$meta.id"
    label 'process_low'

    errorStrategy 'finish'

    // issue with chia, will need work around
    // container "${ workflow.containerEngine == 'singularity' ? 'file:///research/project/shared/benoukraf_lab/.singularity_cache/metascope.sif' : null}"

    cpus 2
    memory '60 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def aligner   = params.ont_long_reads ? "other" : "bowtie"
    """
    cp ${projectDir}/bin/MetaScope_id.R ./
    mkdir bam
    cp -L $bam ./bam/

    singularity run -B ./:\$HOME -C /research/project/shared/benoukraf_lab/.singularity_cache/metascope.sif \\
     Rscript --vanilla MetaScope_id.R \\
        bam/$bam \\
        $aligner
    """
}