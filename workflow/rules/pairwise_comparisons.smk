def get_difftype(x):
	if x in h3k27ac_targets:
		return("h3k27ac")
	else:
		return("binding")


rule create_pairwise_comparisons_rmd:
	input:
		module_pw = "workflow/modules/pairwise_comparison.Rmd",
		r = "workflow/scripts/create_pairwise_comparison.R"
	output:
		rmd = expand(
			os.path.join(
				"analysis", 
				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
				),
			f = "pairwise_comparison.Rmd"
		)
	params:
		git = git_add,
		threads = 4,
	conda: "../envs/rmarkdown.yml"
	retries: git_tries
	threads: 1
	resources:
		runtime = "1m",
		mem_mb = 512
	log: log_path + "/create_rmd/create_{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison_rmd"
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
            {params.threads} \
            {output.rmd} &>> {log}

        ## Add the remainder of the module as literal text
        cat {input.module_pw} >> {output.rmd}

		if [[ {params.git} == "True" ]]; then
			git add {output.rmd}
		fi
		"""

rule compile_pairwise_comparisons_html:
	input:
		annotations = ALL_RDS,
		config = "config/config.yml",
		extrachips = rules.update_extrachips.output,
		here = here_file,
		module_rna = "workflow/modules/rnaseq_pairwise.Rmd",
		results_t1 = lambda wildcards: expand(
			os.path.join(
				"docs", "{{t1}}_{{ref1}}_{{treat1}}_differential_{type}.html"
			),
			type = get_difftype(wildcards.t1)
		),
		results_t2 = lambda wildcards: expand(
			os.path.join(
				"docs", "{{t2}}_{{ref2}}_{{treat2}}_differential_{type}.html"
			),
			type = get_difftype(wildcards.t2)
		),
		rmd = expand(
			os.path.join(
				"analysis", 
				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
				),
			f = "pairwise_comparison.Rmd"
		),
		rmd_config = "config/rmarkdown.yml",
		scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output
	output:
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
		csv = expand(
			os.path.join(
				"output", "pairwise_comparisons", "{{t1}}_{{t2}}",
			"{{t1}}_{{ref1}}_{{treat1}}-{{t2}}_{{ref2}}_{{treat2}}-{f}"
			),
			f = [
				'pairwise_comparison.csv.gz', 'enrichment.csv', 
				'rnaseq_enrichment.csv'
				]
		),
		rds = os.path.join(
			"output", "pairwise_comparisons", "{t1}_{t2}",
			"{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}-all_windows.rds"
		),
		renv = temp(
			os.path.join(
				"output", "envs",
				"{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}-pairwise_comparison.RData"
			)
		)
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = git_tries,
		asset_path = os.path.join(
			"docs", "assets", "{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}"
		)
	conda: "../envs/rmarkdown.yml"
	threads: 4
	resources:
		runtime = "2h",
		mem_mb = 8192
	log: "workflow/logs/pairwise/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison.log"
	shell:
		"""
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}

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
            git add {output.html} {output.fig_path} {output.csv}
			git add {params.asset_path}
        fi
		"""