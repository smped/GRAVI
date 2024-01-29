
args <- commandArgs(TRUE)
target <- args[[1]]
min_prop <- as.numeric(args[[2]])
rmd <- args[[3]]
cat(args, "\n")

library(tidyverse)
library(glue)

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
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
