library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
target <- args[[1]]
threads <- args[[2]]
rmd <- args[[3]]

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
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
