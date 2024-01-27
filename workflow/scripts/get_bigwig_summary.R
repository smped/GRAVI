library(tidyverse)
library(rtracklayer)
library(plyranges)
library(yaml)

stopifnot(library(extraChIPs, logical.return = TRUE))

args <- commandArgs(TRUE)

bw <- args[[1]]
tsv <- args[[2]]
stopifnot(file.exists(bw))

config <- read_yaml(here::here("config", "config.yml"))
sq <- read_rds(
  here::here("output/annotations/seqinfo.rds")
)
blacklist <- import.bed(
  here::here(config$external$blacklist),
  seqinfo = sq
)


gr <- sq %>%
  GRanges() %>%
  setdiff(blacklist)
cov <- import.bw(BigWigFile(bw), which = gr)
cov %>%
  sortSeqlevels() %>%
  split(f = seqnames(.)) %>%
  lapply(function(x){max(x$score)}) %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "seqnames", values_to = "score") %>%
  left_join(as_tibble(sq), by = "seqnames") %>%
  mutate(start = 1) %>%
  dplyr::select(seqnames, start, end = seqlengths, score) %>%
  write_tsv(tsv)
