#!/usr/bin/env Rscript
# do MetaScope aligning against target indices
#
library("MetaScope")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least three argument must be supplied path/to/bam aligner", call.=FALSE)
}

if (file.exists(args[1])) {
	bam = args[1]
}

type = tools::file_ext(args[1])

# Align against target index
metascope_id(input_file = bam, input_type = type, aligner = args[2], num_species_plot = 0)
