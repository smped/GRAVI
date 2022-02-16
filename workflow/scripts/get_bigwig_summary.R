library(tidyverse)
library(rtracklayer)
library(plyranges)

args <- commandArgs(TRUE)

bw <- args[[1]]
tsv <- args[[2]]
stopifnot(file.exists(bw))

sq <- read_rds(
  here::here("output/annotations/seqinfo.rds")
)
blacklist <- import.bed(
  here::here("output/annotations/blacklist.bed.gz"), genome = sq
)
gr <- sq %>%
  GRanges() %>%
  setdiff(blacklist)
cov <- import.bw(BigWigFile(bw), which = gr)
cov %>%
  sortSeqlevels() %>%
  split(f = seqnames(.)) %>%
  lapply(function(x) mutate(range(x), max = max(x$score))) %>%
  GRangesList() %>%
  unlist() %>%
  as.data.frame() %>%
  write_tsv(tsv)
