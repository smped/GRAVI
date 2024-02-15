# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
    conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
    if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
    prev_paths <- .libPaths()
    paths_to_set <- unique(c(conda_lib_path, prev_paths))
    .libPaths(paths_to_set)
}
## A function for printing input
cat_list <- function(x, slot = NULL, sep = "\n\t"){
    nm <- setdiff(names(x), "")
    invisible(
        lapply(
            nm,
            \(i) cat("Received", slot, i, sep, paste0( x[[i]], "\n\t"), "\n")
        )
    )
}

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")

## For testing
# all_input <- list(
#     bam = "../GRAVI_testing/data/bam/SRR8315192.bam",
#     sq = "../GRAVI_testing/output/annotations/seqinfo.rds"
# )
# all_output <- list(
#     bed = "../GRAVI_testing/output/annotations/SRR8315192_greylist.bed.gz"
# )

cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...\n")
library(readr)
library(GreyListChIP)
library(plyranges)
sq <- read_rds(all_input$sq)

## This is set for only a single input file
cat("Initialising new GreyList object\n")
gl <- new("GreyList", karyotype = sq)
cat("Counting reads in ", all_input$bam, "...\n")
gl <- countReads(gl, all_input$bam)
cat("Calculating thresholds...\n")
set.seed(100)
gl <- calcThreshold(gl)
cat("Making greylist...\n")
gl <- makeGreyList(gl)
cat("Writing to ", all_output$bed, "\n")
write_bed(slot(gl, "regions"), all_output$bed)
cat("done\n")
