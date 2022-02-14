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
	editor_options:
	  chunk_output_type: console
	---

	```{r set-knitr-opts, echo=FALSE, child = '../workflow/modules/setup_chunk.Rmd'}
	```

	```{r set-vals}
	target <- \"{{target}}\"
	threads <- {{threads}}
	```

	```{r build-from-module, echo = TRUE, child = '../workflow/modules/macs2_summary.Rmd'}
	```

	",
	.open = "{{",
	.close = "}}"
) %>%
	write_lines(rmd)
