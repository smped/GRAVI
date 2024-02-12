# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}

library(readr)
library(GreyListChIP)

# args <- commandArgs(TRUE)
# bam <- args[[1]]
# sq <- read_rds(args[[2]])
# bed <- args[[3]]
bam <- snakemake@input[["bam"]]
sq <- read_rds(snakemake@input[["seqinfo"]])
bed <- snakemake@output[["bed"]]

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
