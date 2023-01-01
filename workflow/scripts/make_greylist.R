library(readr)
library(GenomeInfoDb)
library(GreyListChIP)
library(Rsamtools)

args <- commandArgs(TRUE)
bam <- BamFile(args[[1]])
sq <- read_rds(args[[2]])
bed <- args[[3]]

## This is set for only a single input file
gl <- new("GreyList", karyotype = sq)
gl <- countReads(gl, bam)
set.seed(100)
gl <- calcThreshold(gl)
gl <- makeGreyList(gl)

export(gl, bed)
