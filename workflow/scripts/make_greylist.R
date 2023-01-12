library(readr)
library(GreyListChIP)

args <- commandArgs(TRUE)
bam <- args[[1]]
sq <- read_rds(args[[2]])
bed <- args[[3]]

## This is set for only a single input file
gl <- new("GreyList", karyotype = sq)
cat("Counting reads in ", bam, "...\n")
gl <- countReads(gl, bam)
cat("Calculating thresholds...\n")
set.seed(100)
gl <- calcThreshold(gl)
cat("Making greylist...\n")
gl <- makeGreyList(gl)
cat("Writing to ", bed, "\n")
export(gl, bed)
cat("done\n")
