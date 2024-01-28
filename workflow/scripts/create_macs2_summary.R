library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
target <- args[[1]]
threads <- args[[2]]
min_prop <- args[[3]]
rmd <- args[[4]]

glue(
	"
	---
	title: '{{target}}: MACS2 Summary'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	    target: \"{{target}}\"
	    threads: {{threads}}
	    min_prop: {{min_prop}}
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
