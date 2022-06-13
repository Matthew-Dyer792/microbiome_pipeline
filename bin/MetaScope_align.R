#!/usr/bin/env Rscript
# do MetaScope aligning against target indices
#
library("MetaScope")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<5) {
  stop("At least three argument must be supplied path/to/fastq path/to/fastq_2 index_prefix output_prefix cpus", call.=FALSE)
}

if (file.exists(args[1])) {
    fastq = args[1]
}

if (args[2]=="null") {
    fastq_2 = NULL
} else {
	if (file.exists(args[2])) {
		fastq_2 = args[2]
	}
}

# Align against target index
align_target_bowtie(
	read1 = fastq,
	read2 = fastq_2,
    lib_dir = "./",
    libs = args[3],
    align_dir = "./",
    align_file = args[4],
    threads = as.numeric(args[5]),
    overwrite = TRUE
)
