
args <- commandArgs(TRUE)
target <- args[[1]]
cores <- as.numeric(args[[2]])
min_prop <- as.numeric(args[[3]])
rmd <- args[[4]]
## I can't understand why args[[2]] is always 1
## The rules show the correct number
cat(args)

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
	    threads: {{cores}}
	    min_prop: {{min_prop}}
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
