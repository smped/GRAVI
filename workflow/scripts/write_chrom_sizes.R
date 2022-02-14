library(tidyverse)
library(GenomicRanges)

args <- commandArgs(TRUE)
infile <- args[[1]]
outfile <- args[[2]]
message("Chrom sizes will be obtained from ", infile)
message("Chrom sizes will be written to ", outfile)

here::here(infile) %>%
  read_tsv(
    col_names = c("chromosome", "start", "end"),
    col_types = "cii--"
  ) %>%
  makeGRangesFromDataFrame(starts.in.df.are.0based = TRUE) %>%
  range() %>%
  as.data.frame() %>%
  dplyr::select(seqnames, end) %>%
  write_tsv(
    here::here(outfile),
    col_names = FALSE
  )
