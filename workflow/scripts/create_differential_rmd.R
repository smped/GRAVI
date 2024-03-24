#' Create the basic header for the differential signal files.
#'
#' Really just need the wildcards to determine all values
#'
#' Also the output will contain the rmd path
#'
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
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

threads <- slot(snakemake, "threads")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_wildcards <- slot(snakemake, "wildcards")
all_params <- slot(snakemake, "params")
cat_list(all_input, "input", sep = ":")
cat_list(all_output, "output", sep = ":")
cat_list(all_wildcards, "wildcards", sep = ":")
cat_list(all_params, "params")

## Solidify file paths
# all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages")
library(tidyverse)
library(glue)

cat_time("Setting main arguments")
target <- all_wildcards$target
ref <- all_wildcards$ref
treat <- all_wildcards$treat
diff_sig_params <- all_params$diff_sig_params
diff_sig_params <- diff_sig_params %>%
    .[vapply(., length, integer(1)) > 0] %>%
    lapply(\(x) ifelse(is.character(x), paste0("\"", x, "\""), x))

cat_time("Writing file header")
glue_params <- names(diff_sig_params) %>%
    lapply(
        \(x) glue("    {x}: {diff_sig_params[[x]]}")
    ) %>%
    .[vapply(., length, integer(1)) > 0] %>%
    c(glue("---\n\n")) %>%
    glue_collapse(sep = "\n")
glue(
    "
	---
	title: '{target} Differential Signal: {treat} Vs. {ref}'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	    counts: \"{all_input$counts}\"
	    treat_levels: [\"{ref}\", \"{treat}\"]"
) %>%
    c(glue_params) %>%
    glue_collapse(sep = "\n") %>%
    write_lines(all_output$rmd)


# glue(
# 	"
# 	---
# 	title: '{{target}} Differential Signal: {{treat}} Vs. {{ref}}'
# 	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
# 	bibliography: references.bib
# 	link-citations: true
# 	params:
# 	    target: \"{{target}}\"
# 	    treat_levels: [\"{{ref}}\", \"{{treat}}\"]
# 	    alpha: \"{{diff_sig_params$alpha}}\"
# 	    fc: \"{{diff_sig_params$fc}}\"
# 	---
#
# 	",
# 	.open = "{{",
# 	.close = "}}"
# ) %>%
# 	write_lines(all_output$rmd)

cat_time("Adding main body of file")
## Now add the rest of the module
file.append(all_output$rmd, all_input$module)
cat_time("Done")
