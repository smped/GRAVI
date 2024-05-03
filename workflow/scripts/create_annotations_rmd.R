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
all_params <- slot(snakemake, "params")

cat_list(all_input, "input:")
cat_list(all_output, "output:")
cat_list(all_params, "params:", "-")

## Solidify file paths
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(glue)
library(yaml)


#### Prepare the RMD ####
cat_time("Forming YAML header")
ln <- glue(
  "
	---
	title: 'Description of Annotations'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	  chrom_sizes: \"{{all_input$chrom_sizes}}\"
	  colours: \"{{all_params$colours}}\"
	  features: \"{{all_input$features}}\"
	  gene_regions: \"{{all_input$gene_regions}}\"
	  gsea_dir: \"{{all_input$gsea_dir}}\"
	  gsea_sig: \"{{all_input$gsea_sig}}\"
	  gtf_exon: \"{{all_input$gtf_exon}}\"
	  gtf_gene: \"{{all_input$gtf_gene}}\"
	  gtf_transcript: \"{{all_input$gtf_transcript}}\"
	  greylist: \"{{all_input$greylist}}\"
	  hic: \"{{all_input$hic}}\"
	  motif_list: \"{{all_input$motifs}}\"
	  motif_uri: \"{{all_input$motif_uri}}\"
	  rna: \"{{all_input$rna}}\"
	  seqinfo: \"{{all_input$seqinfo}}\"
	  trans_models: \"{{all_input$trans_models}}\"
	  tss: \"{{all_input$tss}}\"
	---

	",
  .open = "{{",
  .close = "}}"
)
cat_time("YAML header:")
cat(ln)
readr::write_lines(ln, all_output$rmd)


cat_time("Written YAML header; Appending Module")
file.append(all_output$rmd, all_input$module)

cat_time("Done")

