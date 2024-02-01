def get_greylist_from_target(wildcards):
    ind = df.target == wildcards.target
    return expand(
        os.path.join(annotation_path, "{file}_greylist.bed"),
        file = set(df[ind]['input'])
    )

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
		blacklist = blacklist,
		r = "workflow/scripts/create_macs2_summary.R"
	output:
		rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd")
	params:
		min_prop = lambda wildcards: macs2_qc_param[wildcards.target]['min_prop_reps'],
	conda: "../envs/rmarkdown.yml"
	threads: 1
	log: log_path + "/create_rmd/create_{target}_macs2_summary.log"
	resources:
		mem_mb = 1024,
		runtime = "2m",
	shell:
		"""
		## Create the generic markdown
		Rscript --vanilla \
			{input.r} \
			{wildcards.target} \
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
		config = "config/config.yml",
		cors = os.path.join(
			macs2_path, "{target}", "{target}_cross_correlations.tsv"
		),
		greylist = get_greylist_from_target,
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
			workflow.cores
		)
	resources:
		mem_mb = 8192
	log: log_path + "/macs2_summmary/compile_{target}_macs2_summary.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""
