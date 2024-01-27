library(tidyverse)
library(glue)
library(yaml)

args <- commandArgs(TRUE)
fl <- here::here(args[[1]])

rmd <- read_yaml(here::here("config", "rmarkdown.yml"))

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

glue(
  "```{r setup, echo = FALSE}",
  conda_pre,
  glue_collapse(out, sep = "\n  "),
  "```",
  .open = "{{", .close = "}}",
  .sep = "\n"
) %>%
  write_lines(fl, append = FALSE)
