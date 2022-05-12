rule pairwise_comparisons:
	input:
		annotations = ALL_RDS,
		config = "config/config.yml",
		here = here_file,
		module_pw = "workflow/modules/pairwise_comparison.Rmd",
		module_rna = "workflow/modules/rnaseq_pairwise.Rmd",
		pkgs = rules.install_packages.output,
		r = "workflow/scripts/create_pairwise_comparison.R",
		results_t1 = "docs/{t1}_{ref1}_{treat1}_differential_binding.html",
		results_t2 = "docs/{t2}_{ref2}_{treat2}_differential_binding.html",
		rmd_config = "config/rmarkdown.yml",
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output
	output:
		rmd = expand(
			os.path.join(
				"analysis", 
				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
				),
			f = "pairwise_comparison.Rmd"
		),
		html = expand(
			os.path.join(
				"docs",
				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
			),
			f = "pairwise_comparison.html"
		),
		fig_path = directory(
			os.path.join(
				"docs",
				"{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison_files"
			)
		),
		csv = os.path.join(
			"output", "{t1}_{t2}",
			"{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison.csv.gz"
		)
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	threads: 4
	log: "workflow/logs/pairwise/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison.log"
	shell:
		"""
		## Create the generic markdown header
        Rscript --vanilla \
            {input.r} \
            {wildcards.t1} \
            {wildcards.ref1} \
            {wildcards.treat1} \
			{wildcards.t2} \
            {wildcards.ref2} \
            {wildcards.treat2} \
            {threads} \
            {output.rmd} &>> {log}

        ## Add the remainder of the module as literal text
        cat {input.module_pw} >> {output.rmd}

        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}

        if [[ {params.git} == "True" ]]; then
            TRIES={params.tries}
            while [[ -f .git/index.lock ]]
            do
                if [[ "$TRIES" == 0 ]]; then
                    echo "ERROR: Timeout while waiting for removal of git index lock" &>> {log}
                    exit 1
                fi
                sleep {params.interval}
                ((TRIES--))
            done
            git add {output.rmd}
            git add {output.html}
            git add {output.fig_path}
			git add {output.csv}
        fi
		"""