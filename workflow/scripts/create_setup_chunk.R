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
config <- slot(snakemake, "config")

cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...")
library(tidyverse)
library(glue)
library(yaml)
cat("done\n")

rmd <- read_yaml(all_input$yml)

conda_pre <- "conda_pre <- system2('echo', '$CONDA_PREFIX', stdout = TRUE)
if (conda_pre != \"\") {
  conda_lib_path <- file.path(conda_pre, 'lib', 'R', 'library')
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}"

out <- deparse(rmd$knitr_opts, width.cutoff = 50L) %>%
  str_replace_all("^list\\(", "knitr::opts_chunk$set(\n  ") %>%
  str_trim() %>%
  str_replace_all("\\)$", "\n)")

cat("Writing chunk to", all_output$rmd, "\n")
glue(
  "```{r setup, echo = FALSE}",
  conda_pre,
  glue_collapse(out, sep = "\n  "),
  "```",
  .open = "{{", .close = "}}",
  .sep = "\n"
) %>%
  write_lines(all_output$rmd, append = FALSE)

