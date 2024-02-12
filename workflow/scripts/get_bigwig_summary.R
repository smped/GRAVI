# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}

library(tidyverse)
library(rtracklayer)
library(plyranges)
library(extraChIPs)

bw <- snakemake@input[["bw"]]
tsv <- snakemake@output
sq <- read_rds(snakemake@input[["sq"]])
blacklist <- importPeaks(
  snakemake@input[["blacklist"]], type = "bed", seqinfo = sq
) %>% 
  unlist()


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
