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
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

all_wildcards <- slot(snakemake, "wildcards")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")

cat_list(all_input, "input")
cat_list(all_output, "output")
cat_list(all_wildcards, "wildcards:", "-")
cat_list(all_params, "params:", "-")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(glue)

ln <- glue(
	"
	---
	title: '{{all_wildcards$target}}: MACS2 Summary'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	  annotation_path: \"{{all_params$annotation_path}}\"
	  macs2_path: \"{{all_params$macs2_path}}\"
	  peak_path: \"{{all_params$peak_path}}\"
	  target: \"{{all_wildcards$target}}\"
	---

	",
	.open = "{{",
	.close = "}}"
)
readr::write_lines(ln, all_output$rmd)

cat_time("Written YAML header\nAppending Module\n")

## Now add the rest of the module
file.append(all_output$rmd, all_input$module)

cat_time("Done")
