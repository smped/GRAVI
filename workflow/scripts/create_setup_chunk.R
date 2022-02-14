library(tidyverse)
library(glue)
library(yaml)

args <- commandArgs(TRUE)
fl <- here::here(args[[1]])

rmd <- read_yaml(here::here("config", "rmarkdown.yml"))

out <- deparse(rmd$knitr_opts, width.cutoff = 50L)%>%
  str_replace_all("^list\\(", "knitr::opts_chunk$set(\n  ") %>%
  str_trim() %>%
  str_replace_all("\\)$", "\n)")

glue(
  "```{r setup, echo = FALSE}",
  glue_collapse(out, sep = "\n  "),
  "```",
  .open = "{{", .close = "}}",
  .sep = "\n"
) %>%
  write_lines(fl, append = FALSE)
