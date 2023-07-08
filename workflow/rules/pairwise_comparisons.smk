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
		threads = 4,
	conda: "../envs/rmarkdown.yml"
	threads: 1
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
		"""

rule compile_pairwise_comparisons_html:
	input:
		annotations = ALL_RDS,
		blacklist = blacklist,
		config = "config/config.yml",
		extrachips = rules.update_extrachips.output,
		here = here_file,
		module_rna = "workflow/modules/rnaseq_pairwise.Rmd",
		results_t1 = "docs/{t1}_{ref1}_{treat1}_differential_binding.html",
		results_t2 = "docs/{t2}_{ref2}_{treat2}_differential_binding.html",
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
		asset_path = os.path.join(
			"docs", "assets", "{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}"
		)
	conda: "../envs/rmarkdown.yml"
	threads: 4
	log: "workflow/logs/pairwise/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison.log"
	resources:
		mem_mb = 32768,
		runtime = "2h"
	shell:
		"""
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""