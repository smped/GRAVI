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

library(tidyverse)
library(yaml)
library(glue)
library(magrittr)

args <- commandArgs(TRUE)
fl <- here::here(args[[1]])

config <- read_yaml(here::here("config", "config.yml"))
rmd <- read_yaml(here::here("config", "rmarkdown.yml"))
samples <- read_tsv(here::here(config$samples$file))
# tf_targets <- unique(dplyr::filter(samples, str_to_lower(type) == "tf")$target)
# h3k27ac_targets <- unique(dplyr::filter(samples, str_to_lower(type) == "h3k27ac")$target)
# atac_targets <- unique(dplyr::filter(samples, str_to_lower(type) == "atac")$target)
# all_targets <- c(tf_targets, h3k27ac_targets)
all_targets <- unique(samples$target)
treats <- unique(samples$treat)

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


# ## Sort out the H3K27ac comparisons
# h3k27ac_comparisons <- lapply(
#   config$comparisons$contrasts,
#   function(x) {
#     dplyr::filter(
#       samples, treat %in% x, target %in% h3k27ac_targets
#     ) %>%
#       mutate(treat = factor(treat, levels = x)) %>%
#       distinct(target, treat) %>%
#       arrange(target, treat) %>%
#       group_by(target) %>%
#       summarise(comparison = paste(treat, collapse = "_")) %>%
#       dplyr::filter(comparison == paste(x, collapse = "_")) %>%
#       unite(rmd, everything(), sep ="_", remove = FALSE)
#   }
# ) %>%
#   bind_rows() %>%
#   split(.$target) %>%
#   setNames(c())

# ## Differential H3K27ac YAML
diff_h3k27ac_yaml <- NULL
# if (length(h3k27ac_comparisons)) {
#   diff_h3k27ac_yaml <- list(
#     text = "Differential H3K27ac",
#     menu = h3k27ac_comparisons %>%
#       lapply(
#         function(x){
#           list(
#             text = unique(x$target),
#             menu = lapply(
#               split(x, f = seq_len(nrow(x))),
#               function(y) {
#                 list(
#                   text = str_replace_all(y$comparison, "(.+)_(.+)", "\\2 Vs. \\1"),
#                   href = paste0(y$rmd, "_differential_h3k27ac.html")
#                 )
#               }
#             ) %>%
#               setNames(NULL)
#           )
#         }
#       )
#   )
# }


## Sort out the pairwise comparisons
## This currently automatically finds every possible combination and compares
## them. Alternatives could be manually specifying or manually excluding...
## H3K27ac comparisons are treated identically as TF comparisons at this point
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
      all_targets,
      function(x) {
        list(
          text = x,
          href = glue("{x}_macs2_summary.html")
        )
      }
    )
  ),

  ## Differential TF Signal
  diff_signal_yaml,

  ## Differential H3K27ac Signal
  diff_h3k27ac_yaml,

  ## Pairwise Comparisons
  pairs_yaml

)
site_yaml$navbar$left <- site_yaml$navbar$left[
  vapply(site_yaml$navbar$left, length, integer(1)) > 0
]
other_nav <- setdiff(names(site_yaml$navbar), c("title", "left"))
site_yaml$navbar <- site_yaml$navbar[c("title", "left", other_nav)]

write_yaml(site_yaml, fl)




