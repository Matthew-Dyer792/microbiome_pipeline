#!/usr/bin/env Rscript
# do MetaScope id and reassignment
#
library(MetaScope)
library(magrittr)

NCBI_key <- "522cc5cd1ca59343ba9f55282d9ecfe6c009"
options("ENTREZ_KEY" = NCBI_key)

# setup a variable to accept CLI args
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("two arguments must be supplied path/to/bam aligner", call. = FALSE)
}

if (file.exists(args[1])) {
  bam <- args[1]
}

type <- tools::file_ext(args[1])

# Align against target index
metascope_id(
  input_file = bam,
  input_type = type,
  aligner = args[2],
  num_species_plot = 0
)
