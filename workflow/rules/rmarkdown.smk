rule create_site_yaml:
	input:
		chk = ALL_CHECKS,
		config_yaml = "config/config.yml",
		rmd_yaml = "config/rmarkdown.yml",
		r = "workflow/scripts/create_site_yaml.R",
	output: os.path.join(rmd_path, "_site.yml")
	conda: "../envs/rmarkdown.yml"
	threads: 1
	resources:
		mem_mb = 1024,
		runtime = "5m",
	log: log_path + "/rmarkdown/create_site_yaml.log"
	shell:
		"""
		Rscript --vanilla {input.r} {output} &>> {log}
		"""


rule create_setup_chunk:
	input:
		chk = ALL_CHECKS,
		config = "config/rmarkdown.yml",
		r = "workflow/scripts/create_setup_chunk.R"
	output:
		rmd = "analysis/setup_chunk.Rmd"
	conda: "../envs/rmarkdown.yml"
	threads: 1
	resources:
		mem_mb = 1024,
		runtime = "5m",
	log: log_path + "/rmarkdown/create_setup_chunk.log"
	shell:
		"""
		Rscript --vanilla {input.r} {output.rmd} &>> {log}
		"""

rule create_index_rmd:
	input:
		chk = ALL_CHECKS,
		rmd = os.path.join("workflow", "modules", "index.Rmd"),
	output:
		os.path.join(rmd_path, "index.Rmd")
	threads: 1
	resources:
		mem_mb = 512,
		runtime = "2m",
	shell:
		"""
		cat {input.rmd} > {output}
		"""

rule compile_index_html:
	input:
		html = HTML_OUT,
		rmd = os.path.join(rmd_path, "index.Rmd"),
		setup = rules.create_setup_chunk.output,
		site_yaml = rules.create_site_yaml.output,
		rulegraph = 'workflow/rules/rulegraph.dot'
	output:
		html = "docs/index.html"
	conda: "../envs/rmarkdown.yml"
	threads: 1
	resources:
		mem_mb = 1024,
		runtime = "5m",	
	log: log_path + "/rmarkdown/compile_index_html.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""

rule create_macs2_summary_rmd:
	input:
		chk = ALL_CHECKS,
		module = "workflow/modules/macs2_summary.Rmd",
	output:
		rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd")
	params:
		min_prop = lambda wildcards: macs2_qc_param[wildcards.target]['min_prop_reps'],
		outlier_thresh = lambda wildcards: macs2_qc_param[wildcards.target]['outlier_threshold'],
		macs2_fdr = lambda wildcards: macs2_param[wildcards.target]['fdr'],
	conda: "../envs/rmarkdown.yml"
	threads: 1
	log: log_path + "/create_rmd/create_{target}_macs2_summary.log"
	resources:
		mem_mb = 1024,
		runtime = "2m",
	script:
		"../workflow/scripts/create_macs2_summary.R"


rule compile_macs2_summary_html:
	input:
		annotations = ALL_RDS,
		bw = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}",
				"{{target}}_{treat}_merged_treat_pileup.bw"
			),
			treat = set(df[df.target == wildcards.target]['treat'])
		),
		cors = os.path.join(
			macs2_path, "{target}", "{target}_cross_correlations.tsv"
		),
		consensus_peaks = os.path.join(
			macs2_path, "{target}", "{target}_consensus_peaks.bed.gz"
		),
		greylist = expand(
			os.path.join(annotation_path, "{f}_greylist.bed.gz"),
			f = set(df['input'])
		),
		rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd"),
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output
	output:
		html = "docs/{target}_macs2_summary.html",
		fig_path = directory(
			os.path.join("docs", "{target}_macs2_summary_files", "figure-html")
		),
		renv = temp(
			os.path.join("output", "envs", "{target}_macs2_summary.RData")
		),
	conda: "../envs/rmarkdown.yml"
	threads:
		lambda wildcards: min(
			len(df[df['target'] == wildcards.target]),
			workflow.cores
		)
	resources:
		mem_mb = 8192
	log: log_path + "/macs2_summmary/compile_{target}_macs2_summary.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""
