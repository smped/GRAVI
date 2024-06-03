#' Handle any conda weirdness
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
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%m-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}


log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input:")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(GenomicInteractions)
samples <- here::here(config$samples$file) %>%
  read_tsv()

sq <- all_input$seqinfo %>% read_rds()

#### HiC ####
### Still Needs Checks

cat_time("Preparing HiC Interactions")
hic <- GInteractions()
hic_path <- here::here(config$external$hic)
if (length(hic_path) > 0) {
  if (file.exists(hic_path)) {
    hic <- makeGenomicInteractionsFromFile(hic_path, type = "bedpe")
  }
  cat_time("Parsed", length(hic), "HiC Interactions")
}
stopifnot(is(hic, "GInteractions"))
## Seqinfo objects can be really difficult here. Separate, then reform
keep_int <- (
  seqnames(anchorOne(hic)) %in% seqnames(sq) &
  seqnames(anchorTwo(hic)) %in% seqnames(sq)
)
hic <- hic[keep_int]
hic <- sortSeqlevels(hic)
a1 <- anchorOne(hic) %>% 
  keepStandardChromosomes()
seqlevels(a1) <- seqlevels(sq)
seqinfo(a1) <- sq
a2 <- anchorTwo(hic) %>% 
  keepStandardChromosomes()
seqlevels(a2) <- seqlevels(sq)
seqinfo(a2) <- sq  
hic <- GenomicInteractions(a1, a2)
cat_time("Exporting to", all_output$hic)
write_rds(sort(hic), all_output$hic, compress = "gz")
cat_time("done\n")



