library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
target <- args[[1]]
ref <- args[[2]]
treat <- args[[3]]
threads <- args[[4]]
rmd <- args[[5]]

glue(
	"
	---
	title: '{{target}} Differential Binding: {{treat}} Vs. {{ref}}'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('analysis/setup_chunk.Rmd')}
	```

	```{r set-vals}
	target <- \"{{target}}\"
	treat_levels <- c(\"{{ref}}\", \"{{treat}}\")
	threads <- {{threads}}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
