rule create_site_yaml:
	input:
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		config_yaml = "config/config.yml",
		rmd_yaml = "config/rmarkdown.yml",
		r = "workflow/scripts/create_site_yaml.R",
	output: os.path.join(rmd_path, "_site.yml")
	conda: "../envs/rmarkdown.yml"
	retries: 100
	threads: 1
	log: log_path + "/rmarkdown/create_site_yaml.log"
	shell:
		"""
		Rscript --vanilla {input.r} {output} &>> {log}
		"""


rule create_setup_chunk:
	input:
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		config = "config/rmarkdown.yml",
		r = "workflow/scripts/create_setup_chunk.R"
	output:
		rmd = "analysis/setup_chunk.Rmd"
	conda: "../envs/rmarkdown.yml"
	retries: 100
	threads: 1
	log: log_path + "/rmarkdown/create_setup_chunk.log"
	shell:
		"""
		Rscript --vanilla {input.r} {output.rmd} &>> {log}
		"""

rule create_index_rmd:
	input:
		os.path.join("workflow", "modules", "index.Rmd")
	output:
		os.path.join(rmd_path, "index.Rmd")
	retries: 100
	threads: 1
	shell:
		"""
		cat {input} > {output}
		"""

rule compile_index_html:
	input:
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		html = HTML_OUT,
		rmd = os.path.join(rmd_path, "index.Rmd"),
		setup = rules.create_setup_chunk.output,
		site_yaml = rules.create_site_yaml.output,
		rulegraph = 'workflow/rules/rulegraph.dot'
	output:
		html = "docs/index.html"
	conda: "../envs/rmarkdown.yml"
	threads: 1
	log: log_path + "/rmarkdown/compile_index_html.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""

rule create_macs2_summary_rmd:
	input:
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		module = "workflow/modules/macs2_summary.Rmd",
		blacklist = blacklist,
		r = "workflow/scripts/create_macs2_summary.R"
	output:
		rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd")
	params:
		threads = lambda wildcards: min(
			len(df[df['target'] == wildcards.target]),
			max_threads
		),
		min_prop = lambda wildcards: macs2_qc_param[wildcards.target]['min_prop_reps'],
	conda: "../envs/rmarkdown.yml"
	threads: 1
	log: log_path + "/create_rmd/create_{target}_macs2_summary.log"
	shell:
		"""
		## Create the generic markdown
		Rscript --vanilla \
			{input.r} \
			{wildcards.target} \
			{params.threads} \
			{params.min_prop} \
			{output.rmd} &>> {log}

		## Add the module directly as literal code
		cat {input.module} >> {output.rmd}
		"""


rule compile_macs2_summary_html:
	input:
		annotations = ALL_RDS,
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.{suffix}"),
			sample = set(df[df.target == wildcards.target]['sample']),
			suffix = ['bam', 'bam.bai']
		),
		blacklist = blacklist,
		bw = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}",
				"{{target}}_{treat}_merged_treat_pileup.bw"
			),
			treat = set(df[df.target == wildcards.target]['treat'])
		),
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		config = "config/config.yml",
		cors = os.path.join(
			macs2_path, "{target}", "{target}_cross_correlations.tsv"
		),
		indiv_macs2 = lambda wildcards: expand(
			os.path.join(macs2_path, "{sample}", "{sample}_{suffix}"),
			sample = set(df[df.target == wildcards.target]['sample']),
			suffix = ['callpeak.log', 'peaks.narrowPeak']
		),
		merged_macs2 = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", "{{target}}_{treat}_merged_{suffix}"
			),
			treat = set(df[df.target == wildcards.target]['treat']),
			suffix = ['callpeak.log', 'peaks.narrowPeak']
		),
		qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv"),
		rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd"),
		scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
		setup = rules.create_setup_chunk.output,
		venn_script = os.path.join("workflow", "scripts", "plot_venn.py"),
		yaml = rules.create_site_yaml.output
	output:
		html = "docs/{target}_macs2_summary.html",
		fig_path = directory(
			os.path.join("docs", "{target}_macs2_summary_files", "figure-html")
		),
		peaks = expand(
			os.path.join(macs2_path, "{{target}}", "{{target}}_{file}"),
			file = ['union_peaks.bed', 'treatment_peaks.rds']
		),
		renv = temp(
			os.path.join("output", "envs", "{target}_macs2_summary.RData")
		),
		venn = "docs/assets/{target}/{target}_common_peaks." + fig_type
	params:
		asset_path = os.path.join("docs", "assets", "{target}")
	conda: "../envs/rmarkdown.yml"
	threads:
		lambda wildcards: min(
			len(df[df['target'] == wildcards.target]),
			max_threads
		)
	resources:
		mem_mb = 8192
	log: log_path + "/macs2_summmary/compile_{target}_macs2_summary.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""
