/*
========================================================================================
    Functions
========================================================================================
*/

process SETUP_SAMTOOLS {
    tag "setup_samtools"
    label 'setup'

    if (params.enable_conda) {
        conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    } else {
        container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/samtools:1.14--hb421002_0' : null}"
    }

    script:
    """
    samtools --help
    """
}

process SETUP_MINIMAP2{
    tag "setup_minimap2"
    label 'setup'

    if (params.enable_conda) {
        conda (params.enable_conda ? "bioconda::minimap2=2.24 bioconda::samtools=1.14" : null)
    } else {
        container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' : null}"
    }

    script:
    """
    minimap2 -h
    """
}

process SETUP_PICARD {
    tag "setup_picard"
    label 'setup'

    errorStrategy = { task.exitStatus == 1 ? 'ignore' : 'terminate' } // change ignore to finish for testing

    if (params.enable_conda) {
        conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    } else {
        container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' : null}"
    }

    script:
    command = params.enable_conda ? 'picard' : 'java -jar /usr/local/share/picard-3.0.0-1/picard.jar'
    """
    $command SplitSamByNumberOfReads --help
    """
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SETUP {

    // take: takes nothing as it is only to download the stuff

    main:

    //
    // MODULE: download the packages for the following functions
    //
    SETUP_SAMTOOLS ()

    SETUP_MINIMAP2 ()

    SETUP_PICARD ()

    // emit: outputs nothing as it is only to download the stuff

}

/*
========================================================================================
    THE END
========================================================================================
*/