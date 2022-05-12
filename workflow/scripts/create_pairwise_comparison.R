library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
target1 <- args[[1]]
ref1 <- args[[2]]
treat1 <- args[[3]]
target2 <- args[[4]]
ref2 <- args[[5]]
treat2 <- args[[6]]
threads <- args[[7]]
rmd <- args[[8]]

glue(
  "
	---
	title: '{{target1}}: {{treat1}} Vs. {{ref1}} Compared To {{target2}}: {{treat2}} Vs. {{ref2}}'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	---

	```{r set-knitr-opts, echo=FALSE, child = here::here('workflow/modules/setup_chunk.Rmd')}
	```

	```{r set-vals}
	threads <- {{threads}}
	pairs <- list(
	  \"{{target1}}\" = c(\"{{ref1}}\", \"{{treat1}}\"),
	  \"{{target2}}\" = c(\"{{ref2}}\", \"{{treat2}}\")
	)
	```

	",
  .open = "{{",
  .close = "}}"
) %>%
  write_lines(rmd)