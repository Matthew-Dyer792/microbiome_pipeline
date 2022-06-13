#!/usr/bin/env Rscript
# do MetaScope aligning against target indices
#
library("MetaScope")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least three argument must be supplied path/to/bam index_prefix cpus", call.=FALSE)
}

if (file.exists(args[1])) {
    bam = args[1]
}

# Align against target index
filter_host_bowtie(
	reads_bam = bam,
    lib_dir = "./",
    libs = args[2],
    make_bam = TRUE,
    threads = as.numeric(args[3]),
    overwrite = TRUE
)
