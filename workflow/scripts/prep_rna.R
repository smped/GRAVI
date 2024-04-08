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
threads <- slot(snakemake, "threads")

cat_list(all_input, "input:")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(magrittr)
library(glue)
library(yaml)
library(BiocParallel)
params <- read_yaml(all_input$yaml)
msigdb_params <- params$msigdb
enrich_params <- params$enrichment

cat("Running with", threads, "threads")
bpparam <- MulticoreParam(threads)

#### RNA ####
## Set this up to take multiple files, which can be named in the YAML
## If not named, use the basename as the name with a snipped suffix
## This file should also be checked during the check_files step
rna <- list()
rna_files <- config$external$rna
if (!is.null(rna_files)) {
  cat_time("Importing gtf")
  gtf <- read_rds(all_input$gtf_gene)

  cat_time("Importing", length(rna_files), "RNA datasets")
  if (any(grepl("(xls|xslx|zip|gz)$", rna_files))) {
    cat("RNA can only be provided as uncompressed csv or tsv files")
    stop()
  }
  rna <- bplapply(
    rna_files,
    \(x) {

      ln <- readLines(x, 1)
      sep <- ifelse(grepl("\\t", ln), "\t", ",")
      df <- read.table(x, sep = sep, header = TRUE)
      df <- as_tibble(df)

      ## The key columns are 'gene_id', 'logFC', and 'FDR'
      ## These should be checked earlier
      gn_col <- intersect(
        c("gene_id", "Geneid", "geneid", "ensembl_gene_id", "ensembl_id"),
        names(df)
      )[[1]]
	  if (length(gn_col) == 0) {
		cat("Couldn't detect gene ids in", x)
		stop()
	  }

      fc_col <- intersect(
        c("logFC", "logfc", "lfc", "log2FoldChange"), names(df)
      )[[1]]
	  if (length(fc_col) == 0) {
		cat("Couldn't detect logFC column in", x)
		stop()
	  }

      p_col <- intersect(
        c("PValue", "PVal", "P", "p", "p_value", "p_val", "P.Value", "pvalue"),
        names(df)
      )[[1]]
	  if (length(p_col) == 0) {
		cat("Couldn't detect PValue column in", x)
		stop()
	  }	  

      fdr_col <- intersect(
        c("fdr", "FDR", "adjP", "adj_p", "adj.P.Value", "padj"), names(df)
      )[[1]]
	  if (length(fdr_col) == 0) {
		cat("Couldn't detect FDR column in", x)
		stop()
	  }	  

      df <- dplyr::select(
        df, gene_id = !!sym(gn_col), logFC = !!sym(fc_col),
        PValue = !!sym(p_col), FDR = !!sym(fdr_col)
      )

	  shared_ids <- intersect(df$gene_id, gtf$gene_id)
	  if (length(shared_ids) == 0) {
		cat("RNA-Seq gene ids do not match those in the GTF")
		stop()
	  }
	
	  df

    }, 
	BPPARAM = bpparam
  )
  if (is.null(names(rna)))
    names(rna) <- str_remove_all(basename(rna_files), "\\.(tsv|csv).*$")
}
cat_time("Exporting RNA to", all_output$rna)
write_rds(rna, all_output$rna, compress = "gz")
cat_time("Done")


## GSEA for all RNA-datasets. Preparing here will save time in all downstream
## workflows
rna_gsea_dir <- rna_gsea_sig <- c()
if (!is.null(rna_files)) {

  cat_time("Loading msigdb")
  msigdb <- read_rds(all_input$msigdb)

  cat_time("Loading fgsea")
  library(fgsea)

  cat_time("Preparing GSEA results from RNA")
  enrich_params <- params$enrichment
  gs_list <- msigdb %>%
    split(.$gs_name) %>%
    bplapply(pull, "ensembl_gene", BPPARAM = bpparam)
  cat_time("Preparing directional GSEA results")
  rna_gsea_dir <- rna %>%
    lapply(
      \(x) setNames(-sign(x$logFC) * log10(x$PValue), x$gene_id)
    ) %>%
    lapply(
      \(x) fgseaMultilevel(
        gs_list, x, BPPARAM = bpparam, nPermSimple = 1e4,
        minSize = min(msigdb_params$size)
      )
    ) %>%
    lapply(mutate, padj = p.adjust(pval, enrich_params$adj)) %>%
    lapply(
      dplyr::filter, size >= enrich_params$min_sig, !is.na(pval)
    ) %>%
    lapply(dplyr::select, starts_with("p"), NES, size, leadingEdge)
  cat_time("Preparing non-directional GSEA results")
  rna_gsea_sig <- rna %>%
    lapply(
      \(x) setNames(-log10(x$PValue), x$gene_id)
    ) %>%
    lapply(
      \(x) fgseaMultilevel(
        gs_list, x, BPPARAM = bpparam, nPermSimple = 1e4,
        minSize = min(msigdb_params$size)
      )
    ) %>%
    lapply(mutate, padj = p.adjust(pval, enrich_params$adj)) %>%
    lapply(
      dplyr::filter, size >= enrich_params$min_sig, !is.na(pval)
    ) %>%
    lapply(dplyr::select, starts_with("p"), NES, size, leadingEdge)
}
cat_time("Exporting GSEA results to", all_output$gsea_dir)
write_rds(rna_gsea_dir, all_output$gsea_dir, compress = "gz")
cat_time("Exporting GSEA results to", all_output$gsea_sig)
write_rds(rna_gsea_sig, all_output$gsea_sig, compress = "gz")
cat_time("Done")


