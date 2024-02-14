# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

library(tidyverse)
library(glue)

target <- slot(snakemake, "wildcards")[["target"]]
mod <- slot(snakemake, "input")[["module"]]
rmd <- slot(snakemake, "output")$rmd

macs2_fdr <- slot(snakemake, "params")$macs2_fdr
outlier_thresh <- slot(snakemake, "params")$outlier_thresh
min_prop <- slot(snakemake, "params")$min_prop

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
	write_lines(rmd)

cat("Written YAML header\nAppending Module\n")

## Now add the rest of the module
file.append(rmd, mod)

cat("Done")