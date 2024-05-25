#' This script defines
#'
#' - Seqinfo (sq)
#' - Chromosome sizes (chrom_sizes for bedGraphToBigWig)
#' - Transcript Models for plotting with Gviz
#' - GRanges for genes, transcripts & exons, taken directly from the gtf
#' - TSS
#' - Unique gene-centric regions
#' - Motifs
#'
#' Running this as a stand-alone script removes any dependency on config.yml
#' which reduces the number of times it is re-run by snakemake
#'
#'
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

cat_list(all_input, "input:")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(glue)
library(yaml)
library(msigdbr)
params <- read_yaml(all_input$yaml)

cat_time("Reading gene-level gtf")
gtf <- read_rds(all_input$gtf)

#### MSigDB ####
cat_time("Preparaing MSigDB using msigdbr...")
msigdb_params <- params$msigdb
msigdb <- msigdbr(msigdb_params$species) %>%
  dplyr::filter(
    gs_cat %in% msigdb_params$gs_cat | gs_subcat %in% msigdb_params$gs_subcat,
    ensembl_gene %in% gtf$gene$gene_id
  ) %>%
  dplyr::filter(
    dplyr::n() >= min(msigdb_params$size),
    dplyr::n() <= max(msigdb_params$size),
    .by = gs_name
  )
cat_time("Loaded", length(unique(msigdb$gs_name)), "gene-sets")

cat_time("Updating Gene-Set URLs...")
gs_url <- msigdb %>%
  distinct(gs_cat, gs_subcat, gs_name, gs_url, gs_exact_source) %>%
  mutate(
    gs_url = case_when(
      gs_subcat == "CP:REACTOME" ~ str_remove_all(gs_url, "\\|.+"),
      gs_subcat == "CP:KEGG" ~ paste0("https://www.genome.jp/pathway/", gs_exact_source),
      gs_subcat == "CP:WIKIPATHWAYS" ~ paste0(
        "https://www.wikipathways.org/pathways/", gs_exact_source, ".html"
      ),
      gs_url == "" ~ "http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp",
      TRUE ~ gs_url
    ) %>%
      setNames(gs_name)
  ) %>%
  pull(gs_url)
msigdb$gs_url <- gs_url[msigdb$gs_name]

cat_time("Exporting msigdb to", all_output$msigdb)
write_rds(msigdb, all_output$msigdb, compress = "gz")
cat_time("Done")



