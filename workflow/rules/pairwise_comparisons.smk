rule pairwise_comparisons:
	input:
		annotations = ALL_RDS,
		config = "config/config.yml",
		pkgs = rules.install_packages.output,
		r = "workflow/scripts/create_pairwise_comparison.R",
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output,
		rmd_config = "config/rmarkdown.yml",
		module = "workflow/modules/pairwise_comparison.Rmd"
	output:
		rmd,
		html,
		fig_path,
		bed,
		tsv
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	threads:
	log:
	shell: