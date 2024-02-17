#' @title Update the Navigation Bar
#'
#' @description Update the navigation bar in the _site.yml file
#'
#' @details
#' Uses the information provided in the config file to update the navigation
#' bar
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

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...\n")
library(tidyverse)
library(yaml)
library(glue)
library(magrittr)

cat("Loading data...\n")
rmd <- read_yaml(all_input$yml)
samples <- read_tsv(here::here(config$samples$file))
all_targets <- unique(samples$target)
treats <- unique(samples$treat)

cat("Defining comparisons...\n")
## Sort out the TF comparisons
comparisons <- lapply(
  config$comparisons$contrasts,
  function(x) {
    dplyr::filter(
      samples, treat %in% x, target %in% all_targets
    ) %>%
      mutate(treat = factor(treat, levels = x)) %>%
      distinct(target, treat) %>%
      arrange(target, treat) %>%
      group_by(target) %>%
      summarise(comparison = paste(treat, collapse = "_")) %>%
      dplyr::filter(comparison == paste(x, collapse = "_")) %>%
      unite(rmd, everything(), sep ="_", remove = FALSE)
  }
) %>%
  bind_rows() %>%
  split(.$target) %>%
  setNames(c())

## Differential Signal YAML
cat("Preparing DiffSignal section...\n")
diff_signal_yaml <- NULL
if (length(comparisons)) {
  diff_signal_yaml <- list(
    text = "Differential Signal",
    menu =   comparisons %>%
      lapply(
        function(x){
          list(
            text = unique(x$target),
            menu = lapply(
              split(x, f = seq_len(nrow(x))),
              function(y) {
                list(
                  text = str_replace_all(y$comparison, "(.+)_(.+)", "\\2 Vs. \\1"),
                  href = paste0(y$rmd, "_differential_signal.html")
                )
              }
            ) %>%
              setNames(NULL)
          )
        }
      )
  )
}


## Sort out the pairwise comparisons
## This currently automatically finds every possible combination and compares
## them. Alternatives could be manually specifying or manually excluding...
## H3K27ac comparisons are treated identically as TF comparisons at this point
cat("Checking for viable pairwise comparisons\n")
pairs_yaml <- NULL
all_pairs <- lapply(
  config$comparisons$contrasts,
  function(x) {
    cont = paste(x, collapse = "_")
    target_combs <- map_lgl(
      split(samples, samples$target), \(y) all(x %in% y$treat)
    ) %>%
      which() %>%
      names() %>%
      sort
    paste(target_combs, cont, sep = "_")
  }
) %>%
  unlist() %>%
  sort()
if (length(all_pairs) > 1) {
  # This can only proceed if we have more than one possible comparison
  pairs <- all_pairs %>%
    combn(2) %>%
    t() %>%
    set_colnames(c("c1", "c2")) %>%
    as_tibble() %>%
    mutate(
      pairs = paste(
        str_extract(c1, "^[A-Za-z0-9]+"),
        str_extract(c2, "^[A-Za-z0-9]+"),
        sep = "-"
      ),
      html = paste(c1, c2, "pairwise_comparison.html", sep = "_"),
      across(
        all_of(c("c1", "c2")),
        \(x) str_remove_all(x, "^[A-Za-z0-9]+_")
      ),
      across(
        all_of(c("c1", "c2")),
        \(x) str_replace_all(x, pattern = "(.+)_(.+)", replacement = "\\2 Vs. \\1")
      )
    ) %>%
    unite(comps, c1, c2, sep = " / ") %>%
    dplyr::select(pairs, comps, html)

  cat("Preparing pairwise menu element...\n")
  pairs_yaml <- list(
    text = "Pairwise Comparisons",
    menu = pairs %>%
      split(.$pairs) %>%
      setNames(NULL) %>%
      lapply(
        function(x) {
          list(
            text = unique(x$pairs),
            menu = lapply(
              split(x, f = x$comps),
              function(y) {
                list(
                  text = y$comps,
                  href = y$html
                )
              }
            ) %>%
              setNames(NULL)
          )
        }
      )
  )
}

cat("Finalising yaml structure...\n")
site_yaml <- rmd$rmarkdown_site
site_yaml$navbar$left <- list(
  ## This first item shouldn't change
  list(icon = "fa-home", text = "Home", href = "index.html"),
  ## The Annotations
  list(
    text = "Annotations", href = "annotation_description.html"
  ),
  ## MACS2 Results
  list(
    text = "MACS2 Peak Calling",
    menu = lapply(
      all_targets,
      function(x) {
        list(
          text = x, href = glue("{x}_macs2_summary.html")
        )
      }
    )
  ),

  ## Differential TF Signal
  diff_signal_yaml,

  ## Pairwise Comparisons
  pairs_yaml

)
site_yaml$navbar$left <- site_yaml$navbar$left[
  map_int(site_yaml$navbar$left, length) > 0
]
other_nav <- setdiff(names(site_yaml$navbar), c("title", "left"))
site_yaml$navbar <- site_yaml$navbar[c("title", "left", other_nav)]

cat("Writing output\n")
write_yaml(site_yaml, all_output$yml)
cat("Done")



