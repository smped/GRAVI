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

## Time stamped messages
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

## For testing
# config <- yaml::read_yaml("config/config.yml")
# all_input <- list(
#   samples =config$samples$file,
#   script = here::here("workflow", "scripts", "create_site_yaml.R"),
#   yml = "config/rmarkdown.yml"
# )
# all_output <- list(yml = "analysis/_site.yml")

log <- slot(snakemake, "log")[[1]]
cat_time("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(yaml)
library(glue)
library(magrittr)

cat_time("Loading data...\n")
rmd <- read_yaml(all_input$yml)
cat_time("Defining targets")
samples <- read_tsv(all_input$samples)
all_targets <- unique(samples$target)
treat_by_target <- samples %>%
  split(.$target) %>%
  lapply(pull, "treat") %>%
  lapply(unique)

cat_time("Loading diff_sig_param")
diff_sig_param <- jsonlite::fromJSON(
  here::here("config/json/differential_signal_param.json")
)

## Sort out the TF comparisons
comparisons <- diff_signal_yaml <- NULL
if (length(diff_sig_param)) {
  cat_time("Defining comparisons...\n")
  comparisons <- diff_sig_param %>%
    lapply(pluck, "contrasts") %>%
    lapply(
      matrix, byrow = FALSE, ncol = 2, dimnames = list(c(), c("ref", "treat"))
    ) %>%
    lapply(as_tibble) %>%
    bind_rows(.id = "target") %>%
    mutate(
      valid_treat = lapply(target, \(x) treat_by_target[[x]]),
      keep_ref = mapply(\(x, y) x %in% unlist(y), x = ref, y = valid_treat),
      keep_treat = mapply(\(x, y) x %in% unlist(y), x = treat, y = valid_treat),
    ) %>%
    dplyr::filter(if_all(starts_with("keep"))) %>%
    dplyr::select(-starts_with("keep")) %>%
    mutate(
      comparison = glue("{treat} Vs. {ref}"),
      rmd = glue("{target}_{ref}_{treat}")
    ) %>%
    split(.$target) %>%
    unname()

  if (length(comparisons)) {
    diff_signal_yaml <- list(
      text = "Differential Signal",
      menu = comparisons %>%
        lapply(
          function(x){
            list(
              text = unique(x$target),
              menu = lapply(
                split(x, f = seq_len(nrow(x))),
                function(y) {
                  list(
                    text = as.character(y$comparison),
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
}

## Sort out the pairwise comparisons
## This currently automatically finds every possible combination and compares
## them. Alternatives could be manually specifying or manually excluding...
cat_time("Checking for viable pairwise comparisons\n")
tgt_regex <- paste(all_targets, collapse = "|")
trt_regex <- samples$treat %>%
  unique() %>%
  paste(collapse = "|")
pairs_yaml <- NULL
if (length(comparisons) > 1) {
  cat_time("Preparing pairwise YAML section")
  pairs_yaml <- list(
    text = "Pairwise Comparisons",
    menu = comparisons %>%
      bind_rows() %>%
      pull("rmd") %>%
      sort() %>%
      combn(2) %>%
      t() %>%
      set_colnames(c("rmd1", "rmd2")) %>%
      as_tibble() %>%
      mutate(
        tgt1 = str_extract(rmd1, paste0("^(", tgt_regex, ")_")),
        tgt2 = str_extract(rmd2, paste0("^(", tgt_regex, ")_")),
        across(starts_with("tgt"), \(x) str_remove_all(x, "_$")),
        cont1 = rmd1 %>%
          str_remove_all(paste0("^(", paste(all_targets, collapse = "|"), ")_")) %>%
          str_replace_all(paste0("^(", trt_regex, ")_(", trt_regex, ")$"), "\\2 Vs. \\1"),
        cont2 = rmd2 %>%
          str_remove_all(paste0("^(", paste(all_targets, collapse = "|"), ")_")) %>%
          str_replace_all(paste0("^(", trt_regex, ")_(", trt_regex, ")$"), "\\2 Vs. \\1"),
        ref = paste0(rmd1, "-", rmd2, "_pairwise_comparison.html")
      ) %>%
      unite(text, starts_with("cont"), sep = " / ") %>%
      unite(menu, starts_with("tgt"), sep = "-") %>%
      distinct(ref, .keep_all = TRUE) %>%
      split(.$menu) %>%
      lapply(
        \(x) {
          list(
            text = unique(x$menu),
            menu = lapply(
              seq_len(nrow(x)),
              \(i) {
                list(
                  text = x$text[[i]],
                  href = x$ref[[i]]
                )
              }
            )
          )
        }
      ) %>%
      unname()
  )
}

cat_time("Checking for additional modules...")
module_yaml <- NULL
add_modules <- vapply(config$modules, \(x) any(x %in% all_targets), logical(1))
if (any(add_modules)) {
  module_yaml <- add_modules %>%
      which() %>%
      names() %>%
      lapply(
        \(x) {
          list(
            text = str_to_upper(x),
            menu = lapply(
              config$modules[[x]],
              \(i) list(text = i, href = glue("{i}_{x}.html"))
            )
          )
        }
      )
}


cat_time("Finalising yaml structure...\n")
shared <- NULL
if (length(all_targets) > 1) {
  shared <- list(list(text = "All Targets", href = "signal_comparison.html"))
}
site_yaml <- rmd$rmarkdown_site
site_yaml$navbar$left <- list(
  ## This first item shouldn't change
  list(icon = "fa-home", text = "Home", href = "index.html"),
  ## The Annotations
  list(
    text = "Annotations", href = "annotation_description.html"
  ),
  ## MACS2 Results
  ## NFR targets (or maybe ROSE eventually)
  module_yaml[[1]], # Needs to change if ROSE is added

  list(
    text = "Signal Detection",
    menu = lapply(
      all_targets,
      function(x) {
        list(
          text = x, href = glue("{x}_signal_summary.html")
        )
      }
    ) %>%
    c(shared)
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

cat(as.yaml(site_yaml))

cat_time("Writing output\n")
write_yaml(site_yaml, all_output$yml)
cat_time("Done")



