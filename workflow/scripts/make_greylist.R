library(readr)
library(GenomeInfoDb)
library(GreyListChIP)

args <- commandArgs(TRUE)
bam <- args[[1]] ## Needs to be the full filename with .bam suffix
sq_file <- args[[2]]
bed <- args[[3]]

sq <- read_rds(sq_file)

## This is set for only a single input file
ip <- file.path(bam_path, "Input", bam)
gl <- new("GreyList", karyotype = sq)
gl <- countReads(gl, bam)
set.seed(100)
gl <- calcThreshold(gl)
gl <- makeGreyList(gl)

export(gl, bed)
