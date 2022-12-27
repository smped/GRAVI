library(tidyverse)
library(yaml)
library(GreyListChIP)

args <- commandArgs(TRUE)
bam <- args[[1]] ## Needs to be the full filename with .bam suffix

config <- read_yaml(here::here("config", "config.yml"))
bam_path <- here::here(config$paths$bam)


## This is set for only a single input file
ip <- file.path(bam_path, "Input", bam)
gl <- new("GreyList", karyotype = sq)
gl <- countReads(gl)
set.seed(100)
gl <- calcThreshold(gl)
gl <- makeGreyList(gl)

bed <- here::here(
  "output", "annotations",
  str_replace(bam, ".bam", "_greylist.bed")
)
export(gl, bed)
