/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def valid_params = [
//     ena_metadata_fields : ['run_accession', 'experiment_accession', 'library_layout', 'fastq_ftp', 'fastq_md5']
// ]

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowSra.initialise(params, log, valid_params)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { PYCOQC                                } from '../modules/pycoqc/main'
include { MINIMAP2_ALIGN                        } from '../modules/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_FILTER     } from '../modules/minimap2/align/main'
include { SAMTOOLS_MERGE as FILTERED_BAM_MERGE  } from '../modules/samtools/merge/main'
include { SAMTOOLS_FASTQ as BAM_TO_FASTQ        } from '../modules/samtools/fastq/main'
include { SAMTOOLS_VIEW as FILTER_BAM           } from '../modules/samtools/view/main'
include { QNAMES as TARGET_QNAMES               } from '../modules/samtools/qnames/main'
include { QNAMES as FILTER_QNAMES               } from '../modules/samtools/qnames/main'
include { MERGE_QNAMES as MERGE_TARGET_QNAMES   } from '../modules/command_line/merge_qnames/main'
include { MERGE_QNAMES as MERGE_FILTER_QNAMES   } from '../modules/command_line/merge_qnames/main'
include { GENERATE_FILTERED_QNAMES              } from '../modules/command_line/compare_qnames/main'
include { SAMTOOLS_INDEX                        } from '../modules/samtools/index/main'
include { BAM_MARKDUPLICATES_PICARD             } from '../subworkflows/bam_markduplicates_picard/main'
include { BAM_STATS_SAMTOOLS                    } from '../subworkflows/bam_stats_samtools/main'
include { MULTIQC                               } from '../modules/multiqc/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// channel for the target genome index files
Channel
    .fromPath(params.target_index, type: 'file', checkIfExists: true)
    .set{ target_index }

// channel for the filter genome index files
Channel
    .fromPath(params.filter_index, type: 'file', checkIfExists: true)
    .set { filter_index }

// empty channel for the cram fasta samtools view requires
channel
    .of( tuple([id: "blah"], []) )
    .first()
    .set { samtools_fasta }

// empty channel for the fai picard markduplicates requires
channel
    .of( tuple([id: "blah"], []) )
    .first()
    .set { picard_fai }

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow ONT_LONG_READS {

    take:
    fastq_files // channel: val( meta ), path([ fastq ])

    main:

    //
    // Setup the target indexs along with the files to be aligned 
    //

    // channel for composition of fastq reads and target indexed genomes
    fastq_files
        .combine(target_index)
        .map{ it-> 
            tuple([id: "${it[0].id}", target_index_prefix: "${it[0].index_prefix}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            params.single_end ? file("${it[1][0]}") : [file("${it[1][0]}"), file("${it[1][1]}")],
            file("${it[2]}")) }
        .set{ target_align_files }



    // look at adding a section for the duplicating this with kraken/centrifuge



    if (params.sequence_summary) {
        // channel for the ONT sequencing_summary.txt files
        Channel
            Channel
            .fromFilePairs(params.sequence_summary, size: 1, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find any .txt files matching: ${params.sequence_summary}\nNB: Path needs to be enclosed in quotes!\n" }
            .map { it -> tuple([id: "${it[0]}"], it[1]) }
            .set{ summary_files }

        //
        // MODULE: Align the raw fastq files against the given microbe index files
        //
        PYCOQC (summary_files)
    }

    //
    // MODULE: Align the raw fastq files against the given microbe index files
    //
    MINIMAP2_ALIGN (target_align_files, true)

    //
    // MODULE: Generate fastq files from the aligned Bam files in the previous step
    //
    BAM_TO_FASTQ (MINIMAP2_ALIGN.out.bam, false)

    //
    // MODULE: Extract the qnames from the target aligned reads
    //
    TARGET_QNAMES (MINIMAP2_ALIGN.out.bam)

    // refactor id to allow grouping of qname files by sample id
    TARGET_QNAMES.out.qname
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ target_qname }

    //
    // MODULE: Merge the target qnames and then return a file containing only the unique target qnames
    //
    MERGE_TARGET_QNAMES (target_qname)

    // cartesian mulitply the sample fastqs to the indexed filter genomes
    BAM_TO_FASTQ.out.other
        .combine(filter_index)  // does this work for a multiple filter indexs
        .map{ it-> 
            tuple([id: "${it[0].id}", target_index_prefix: "${it[0].index_prefix}", single_end: params.single_end, index_prefix: "${it[2].simpleName}"],
            file("${it[1]}"),
            file("${it[2]}")) }
        .set{ filter_files }

    //
    // MODULE: Align the converted fastq files against the given filter index files (human, etc...)
    //
    MINIMAP2_FILTER (filter_files, true)

    //
    // MODULE: Extract the qnames from the target aligned reads
    //
    FILTER_QNAMES (MINIMAP2_FILTER.out.bam)

    // refactor id to allow grouping of bam files by sample id
    FILTER_QNAMES.out.qname
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ filter_qnames }

    //
    // MODULE: Merge the filter qnames and then return a file containing only the unique filter qnames
    //
    MERGE_FILTER_QNAMES (filter_qnames)

    // refactor id to allow grouping of qname files by sample id
    MERGE_TARGET_QNAMES.out.qname
        .mix( MERGE_FILTER_QNAMES.out.qname )
        .groupTuple()
        .set{ filtered_qname }

    //
    // MODULE: Merge the target & filter qnames together and then return a file containing only the qnames present in both
    //
    GENERATE_FILTERED_QNAMES (filtered_qname)

    // need to remap to set the 2 index of the tuple to a shared ID
    MINIMAP2_ALIGN.out.bam
        .map{ it->
            tuple([id: "${it[0].id}", single_end: params.single_end, index_prefix: "${it[0].index_prefix}"],
            file("${it[1]}"),
            "${it[0].id}") }
        .set{ target_bam }

    // need to remap to set the 2 index of the tuple to a shared ID
    GENERATE_FILTERED_QNAMES.out.qname
        .map{ it->
            tuple([id: "${it[0].id}"],
            file("${it[1]}"),
            "${it[0].id}") }
        .set{ filtered_qnames }

    // cartesian mulitply the aligned bams to the filtered qnames
    target_bam
        .combine( filtered_qnames, by: 2 ) // 2 is the index of the shared ID
        .map{ it-> 
            tuple([id: "${it[1].id}", single_end: params.single_end, index_prefix: "${it[1].index_prefix}"],
            file("${it[2]}"),
            file("${it[4]}")) }
        .set{ to_filter }

    //
    // MODULE: Filter the target bam file to remove all reads that aligned to the filter genomes
    //
    FILTER_BAM (to_filter, samtools_fasta)

    // refactor id to allow grouping of filtered bam files by sample id
    FILTER_BAM.out.bam
        .map{ it -> tuple([id: "${it[0].id}"], it[1]) }
        .groupTuple()
        .set{ to_merge }

    //
    // MODULE: Merge the filtered bam files together
    //
    FILTERED_BAM_MERGE (to_merge)

    //
    // SUBWORKFLOW: Sort, Index, Mark duplicates, and Stats on the merged bam
    //
    BAM_MARKDUPLICATES_PICARD (FILTERED_BAM_MERGE.out.bam, samtools_fasta, picard_fai)

    //
    // MODULE: Remove duplicates
    //
    

    // prepare the channel of target and filter alignments
    MINIMAP2_ALIGN.out.bam
        .mix( MINIMAP2_FILTER.out.bam )
        .set{ bam_to_index }

    // bam_to_stats.view()

    //
    // MODULE: Index all aligned files
    //
    SAMTOOLS_INDEX ( bam_to_index )

    // merge the sorted bam and index files together by id
    bam_to_index
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    // ch_bam_bai.view()

    //
    // SUBWORKFLOW: Produce stats on the target and filter alignments
    //
    BAM_STATS_SAMTOOLS (ch_bam_bai, samtools_fasta)

    // create a channel containing the multiqc files
    BAM_STATS_SAMTOOLS.out.stats
        .mix(BAM_STATS_SAMTOOLS.out.flagstat, BAM_STATS_SAMTOOLS.out.idxstats, BAM_MARKDUPLICATES_PICARD.out.metrics, BAM_MARKDUPLICATES_PICARD.out.stats, BAM_MARKDUPLICATES_PICARD.out.flagstat, BAM_MARKDUPLICATES_PICARD.out.idxstats, PYCOQC.out.json)
        .map { it -> it[1] }
        .collect()
        .set{ mulitqc_files }

    //
    // STEP: Summarize results with multiqc
    //
    MULTIQC (mulitqc_files, params.multiqc_config, params.extra_multiqc_config, params.multiqc_logo)


    emit:
    bam         = FILTERED_BAM_MERGE.out.bam      // channel: [ val(meta), bam   ]

}

/*
========================================================================================
    THE END
========================================================================================
*/