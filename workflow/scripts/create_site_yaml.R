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
pairs <- config$pairwise %>%
  lapply(
    function(x) {
      tibble(
        pairs = paste(names(x), collapse = "-"),
        comps = vapply(
          x,
          function(y) paste(rev(y), collapse = " Vs. "),
          character(1)
        ) %>%
          paste(collapse = " / "),
        rmd = paste(
          vapply(
            names(x),
            function(y) paste(y, paste(x[[y]], collapse = "_"), sep = "_"),
            character(1)),
          collapse = "_") %>%
          paste0("_pairwise_comparison.Rmd")
      )
    }
  ) %>%
  bind_rows()

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

annotation_html <- list(
      list(
        text = "Genome and Transcriptome Annotations",
        href = "annotation_description.html"
      )
)

site_yaml$navbar$left <- list(
  ## This first item shouldn't change
  list(icon = "fa-home", text = "Home", href = "index.html"),
  ## The Annotations
  list(
    text = "Annotations",
    menu = annotation_html
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
            text = x$target,
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
  list(
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
                  href = y$rmd
                )
              }
            ) %>%
              setNames(NULL)
          )
        }
      )
  )

)
other_nav <- setdiff(names(site_yaml$navbar), c("title", "left"))
site_yaml$navbar <- site_yaml$navbar[c("title", "left", other_nav)]

write_yaml(site_yaml, fl)




