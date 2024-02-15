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

all_wildcards <- slot(snakemake, "wildcards")
all_params <- slot(snakemake, "params")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")

cat_list(all_wildcards, "wildcards", ":")
cat_list(all_params, "params")
cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...\n")
library(tidyverse)
library(glue)

target <- all_wildcards$target
macs2_fdr <- all_params$macs2_fdr
outlier_thresh <- all_params$outlier_thresh
min_prop <- all_params$min_prop

glue(
	"
	---
	title: '{{target}}: MACS2 Summary'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	    target: \"{{target}}\"
	    min_prop: {{min_prop}}
	    outlier_thresh: {{outlier_thresh}}
	    macs2_fdr: {{macs2_fdr}}
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(all_output$rmd)

cat("Written YAML header\nAppending Module\n")

## Now add the rest of the module
file.append(all_output$rmd, all_input$module)

cat("Done")
