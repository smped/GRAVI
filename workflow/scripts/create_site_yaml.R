#' @title Update the Navigation Bar
#'
#' @description Update the navigation bar in the _site.yml file
#'
#' @details
#' Uses the information provided in the config file to update the navigation
#' bar
library(tidyverse)
library(yaml)
library(glue)
library(magrittr)

args <- commandArgs(TRUE)
fl <- here::here(args[[1]])

config <- read_yaml(here::here("config", "config.yml"))
rmd <- read_yaml(here::here("config", "rmarkdown.yml"))
samples <- read_tsv(here::here(config$samples$file))
targets <- sort(unique(samples$target))
treats <- unique(samples$treat)

## Sort out the comparisons
comparisons <- lapply(
  config$comparisons$contrasts,
  function(x) {
    dplyr::filter(samples, treat %in% x) %>%
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

## Sort out the pairwise comparisons
## This currently automatically finds every possible combination and compares
## them. Alternatives could be manually specifying or manually excluding...

pairs_yaml <- NULL
all_pairs <- lapply(
  config$comparisons$contrasts,
  function(x) {
    cont = paste(x, collapse = "_")
    target_combs <- vapply(
      split(samples, samples$target),
      function(y) all(x %in% y$treat),
      logical(1)
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
        str_remove_all, "^[A-Za-z0-9]+_"
      ),
      across(
        all_of(c("c1", "c2")),
        str_replace_all, pattern = "(.+)_(.+)", replacement = "\\2 Vs. \\1"
      )
    ) %>%
    unite(comps, c1, c2, sep = " / ") %>%
    dplyr::select(pairs, comps, html)
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

site_yaml <- rmd$rmarkdown_site

if (!is.null(site_yaml$output_dir)) {
  out_dir <- file.path(dirname(fl), site_yaml$output_dir)
  message("Checking for directory: ", out_dir)
  if (!dir.exists(out_dir)) {
    message("Creating: ", out_dir)
    stopifnot(dir.create(out_dir))
  }
  message(out_dir, " exists")
}

site_yaml$navbar$left <- list(
  ## This first item shouldn't change
  list(icon = "fa-home", text = "Home", href = "index.html"),
  ## The Annotations
  list(
    text = "Annotations",
    href = "annotation_description.html"
  ),
  ## MACS2 Results
  list(
    text = "MACS2 Peak Calling",
    menu = lapply(
      targets,
      function(x) {
        list(
          text = x,
          href = glue("{x}_macs2_summary.html")
        )
      }
    )
  ),
  ## Differential Binding
  list(
    text = "Differential Binding",
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
                  href = paste0(y$rmd, "_differential_binding.html")
                )
              }
            ) %>%
              setNames(NULL)
          )
        }
      )
  ),

  ## Pairwise Comparisons
  pairs_yaml

)
site_yaml$navbar$left <- site_yaml$navbar$left[
  vapply(site_yaml$navbar$left, length, integer(1)) > 0
]
other_nav <- setdiff(names(site_yaml$navbar), c("title", "left"))
site_yaml$navbar <- site_yaml$navbar[c("title", "left", other_nav)]

write_yaml(site_yaml, fl)




