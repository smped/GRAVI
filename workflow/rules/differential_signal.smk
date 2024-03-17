# rule create_differential_signal_rmd:
# 	input:
# 		chk = ALL_CHECKS,
# 		module = os.path.join("workflow", "modules", "differential_signal.Rmd"),
# 		r = os.path.join("workflow", "scripts", "create_differential_rmd.R")
# 	output: 
# 		rmd = os.path.join(
# 			rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
# 		)
# 	conda: "../envs/rmarkdown.yml"
# 	params:
# 		win_type = "fixed" # Add the rest really soon!!!
# 	log: log_path + "/create_rmd/{target}_{ref}_{treat}_differential_signal.log"
# 	threads: 1
# 	shell:
# 		"""
# 		## Create the generic markdown header
# 		Rscript --vanilla \
# 			{input.r} \
# 			{wildcards.target} \
# 			{wildcards.ref} \
# 			{wildcards.treat} \
# 			{output.rmd} &>> {log}

# 		## Add the remainder of the module as literal text
# 		cat {input.module} >> {output.rmd}
# 		"""

rule count_windows:
	input:
		bam = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.bam"),
			sample = df['sample'][(df['target'] == wildcards.target)]
		),
		bai = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.bam.bai"),
			sample = df['sample'][(df['target'] == wildcards.target)]
		),
		input_bam = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.bam"),
			sample = set(df['input'][(df['target'] == wildcards.target)])
		),
		input_bai = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.bam.bai"),
			sample = set(df['input'][(df['target'] == wildcards.target)])
		),
		blacklist = os.path.join(annotation_path, "blacklist.rds"),
		chk = ALL_CHECKS,
		greylist = lambda wildcards: expand(
			os.path.join(grey_path, "{ip_sample}_greylist.bed.gz"),
			ip_sample = set(df['input'][df['target'] == wildcards.target])
		),
		macs2_logs = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", 
				"{{target}}_{treat_levels}_merged_callpeak.log"
			),
			treat_levels = set(df['treat'][df['target'] == wildcards.target])
		),
		macs2_qc = os.path.join(
			macs2_path, "{target}", "{target}_qc_samples.tsv"
		),
		peaks = lambda wildcards: expand(
			os.path.join(
				peak_path, "{{target}}",
				"{{target}}_{treat}_filtered_peaks.narrowPeak"
			),
			treat = set(df['treat'][df['target'] == wildcards.target])
		),
		script = os.path.join("workflow", "scripts", "make_counts.R"),
		seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
	output:
		rds = os.path.join("data", "counts", "{target}_counts.rds")
	conda: "../envs/rmarkdown.yml"
	log: os.path.join(log_path, "count_windows", "{target}_make_counts.log")
	threads: 8
	params:
		contrasts = lambda wildcards: diff_sig_param[wildcards.target]['contrasts'],
		filter_q = lambda wildcards: diff_sig_param[wildcards.target]['filter_q'],
		win_type = lambda wildcards: diff_sig_param[wildcards.target]['window_type'],
		win_size = lambda wildcards: diff_sig_param[wildcards.target]['window_size'],
		win_step = lambda wildcards: diff_sig_param[wildcards.target]['window_step'],
	resources:
		runtime = "1h",
		mem_mb = 64000,
	script:
		"../scripts/make_counts.R"

# rule compile_differential_signal_html:
# 	input:
# 		annotations = ALL_RDS,
# 		merged_macs2 = lambda wildcards: expand(
# 			os.path.join(
# 				macs2_path, "{{target}}", "{{target}}_{pre}_merged_callpeak.log"
# 			),
# 			pre = [wildcards.ref, wildcards.treat]
# 		),
# 		merged_bw = lambda wildcards: expand(
# 			os.path.join(
# 				macs2_path, "{{target}}", 
# 				"{{target}}_{pre}_merged_treat_pileup.bw"
# 			),
# 			pre = [wildcards.ref, wildcards.treat]
# 		),
# 		peaks = expand(
# 			os.path.join(macs2_path, "{target}", "{target}_union_peaks.bed"),
# 			target = targets
# 		),
# 		module = os.path.join("workflow", "modules", db_method + ".Rmd"),
# 		rmd = os.path.join(
# 			rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
# 		),
# 		samples = os.path.join(
# 			macs2_path, "{target}", "{target}_qc_samples.tsv"
# 		),
# 		setup = rules.create_setup_chunk.output,
# 		site_yaml = rules.create_site_yaml.output,
# 		yml = expand(
# 			os.path.join("config", "{file}.yml"), file = ['config', 'params']
# 		),
# 		rnaseq_mod = os.path.join(
# 			"workflow", "modules", "rnaseq_differential.Rmd"
# 		),
# 		windows = rules.count_windows.output
# 	output:
# 		html = "docs/{target}_{ref}_{treat}_differential_signal.html",
# 		fig_path = directory(
# 			os.path.join(
# 				"docs", "{target}_{ref}_{treat}_differential_signal_files",
# 				"figure-html"
# 			)
# 		),
# 		renv = temp(
# 			os.path.join(
# 				"output", "envs",
# 				"{target}_{ref}_{treat}-differential_signal.RData"
# 			)
# 		),
# 		outs = expand(
# 			os.path.join(
# 				diff_path, "{{target}}", "{{target}}_{{ref}}_{{treat}}-{f}"
# 			),
# 			f = [
# 				'differential_signal.rds', 'down.bed', 'up.bed','differential_signal.csv.gz', 'DE_genes.csv', 'enrichment.csv',
# 				'rnaseq_enrichment.csv'
# 			]
# 		)
# 	retries: 1
# 	conda: "../envs/rmarkdown.yml"
# 	threads:
# 		lambda wildcards: min(
# 			len(
# 				df[
# 					(df['target'] == wildcards.target) &
# 					(
# 						(df['treat'] == wildcards.ref) |
# 						(df['treat'] == wildcards.treat)
# 					)
# 				]
# 			),
# 			workflow.cores
# 		)
# 	resources:
# 		mem_mb = 65536,
# 		runtime = "3h"
# 	log: log_path + "/differential_signal/{target}_{ref}_{treat}_differential_signal.log"
# 	shell:
# 		"""
# 		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
# 		"""
