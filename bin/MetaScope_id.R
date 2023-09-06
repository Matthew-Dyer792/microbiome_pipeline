#!/usr/bin/env Rscript
# do MetaScope id and reassignment
#
library(MetaScope)
library(magrittr)

# setup a variable to accept CLI args
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("two arguments must be supplied path/to/bam aligner", call. = FALSE)
}

NCBI_key <- args[1]
options("ENTREZ_KEY" = NCBI_key)

if (file.exists(args[2])) {
  bam <- args[2]
}

type <- tools::file_ext(bam)

# Align against target index
metascope_id(
  input_file = bam,
  input_type = type,
  aligner = args[3],
  num_species_plot = 0
)

# close by writing out sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
