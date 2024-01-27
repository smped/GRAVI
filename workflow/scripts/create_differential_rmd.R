library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
target <- args[[1]]
ref <- args[[2]]
treat <- args[[3]]
threads <- args[[4]]
type <- args[[5]]
rmd <- args[[6]]

glue(
	"
	---
	title: '{{target}} Differential {{type}}: {{treat}} Vs. {{ref}}'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	    target: \"{{target}}\"
	    treat_levels: [\"{{ref}}\", \"{{treat}}\"]
	    threads: {{threads}}
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
