library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
target1 <- args[[1]]
target2 <- args[[2]]
ref <- args[[3]]
treat <- args[[4]]
threads <- args[[5]]
rmd <- args[[6]]

glue(
  "
	---
	title: '{{treat}} Vs. {{ref}}: Comparison between {{target1}} and {{target2}}'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	---

	```{r set-knitr-opts, echo=FALSE, child = '../workflow/modules/setup_chunk.Rmd'}
	```

	```{r set-vals}
	targets <- c(\"{{target1}}\", \"{{target2}}\")
	treat_levels <- c(\"{{ref}}\", \"{{treat}}\")
	threads <- {{threads}}
	```

	```{r build-from-module, echo = TRUE, child = '../workflow/modules/pairwise_comparison.Rmd'}
	```

	",
  .open = "{{",
  .close = "}}"
) %>%
  write_lines(rmd)
