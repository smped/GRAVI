rule make_greylist:
	input: 
		bam = os.path.join(bam_path, "{ip_sample}.bam"),
		bim = os.path.join(bam_path, "{ip_sample}.bam.bai"),
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		r = os.path.join("workflow", "scripts", "make_greylist.R"),
		seqinfo = os.path.join(annotation_path, "seqinfo.rds")
	output:
		bed = os.path.join(annotation_path, "{ip_sample}_greylist.bed")
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/scripts/{ip_sample}_make_greylist.log"
	threads: 1
	resources:
		mem_mb = 16384
	shell:
		"""
		Rscript --vanilla \
			{input.r} \
			{input.bam} \
			{input.seqinfo} \
			{output.bed} &>> {log}
		"""

rule create_differential_signal_rmd:
	input:
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		db_mod = os.path.join(
			"workflow", "modules", "differential_signal.Rmd"
		),
		r = "workflow/scripts/create_differential_rmd.R"
	output: 
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
		)
	params:
		threads = lambda wildcards: min(
			len(
				df[
					(df['target'] == wildcards.target) &
					(
						(df['treat'] == wildcards.ref) |
						(df['treat'] == wildcards.treat)
					)
				]
			),
			max_threads
		),
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/create_rmd/{target}_{ref}_{treat}_differential_signal.log"
	threads: 1
	shell:
		"""
		## Create the generic markdown header
		Rscript --vanilla \
			{input.r} \
			{wildcards.target} \
			{wildcards.ref} \
			{wildcards.treat} \
			{params.threads} \
			{output.rmd} &>> {log}

		## Add the remainder of the module as literal text
		cat {input.db_mod} >> {output.rmd}
		"""

rule count_windows:
	input:
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.{suffix}"),
			sample = df['sample'][(df['target'] == wildcards.target)],
			suffix = ['bam', 'bam.bai']
		),
		blacklist = blacklist,
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		greylist = lambda wildcards: expand(
			os.path.join(annotation_path, "{ip_sample}_greylist.bed"),
			ip_sample = set(df['input'][df['target'] == wildcards.target])
		),
		peaks = os.path.join(
			macs2_path, "{target}", "{target}_union_peaks.bed"
		),
		samples = os.path.join(
			macs2_path, "{target}", "{target}_qc_samples.tsv"
		),
		r = os.path.join("workflow", "scripts", "make_filtered_counts.R")
	output:
		rds = expand(
			os.path.join(
				"differential_signal", "{{target}}", "{{target}}_{file}.rds"
			),
			file = ['window_counts', 'filtered_counts']
		)
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/diferential_signal/{target}_count_windows.log"
	threads:
		lambda wildcards: min(
			len(df[(df['target'] == wildcards.target)]), max_threads
		)
	resources:
		runtime = "2h",
		mem_mb = 32768
	shell:
		"""
		Rscript --vanilla \
		  {input.r} \
		  {wildcards.target} \
		  {threads} &>> {log}
		"""

rule compile_differential_signal_html:
	input:
		annotations = ALL_RDS,
		chk = expand(
			os.path.join("output", "checks", "{f}.chk"),
			f = ['r-packages', 'here']
		),
		merged_macs2 = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", "{{target}}_{pre}_merged_callpeak.log"
			),
			pre = [wildcards.ref, wildcards.treat]
		),
		merged_bw = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", 
				"{{target}}_{pre}_merged_treat_pileup.bw"
			),
			pre = [wildcards.ref, wildcards.treat]
		),
		peaks = expand(
			os.path.join(macs2_path, "{target}", "{target}_union_peaks.bed"),
			target = targets
		),
		module = os.path.join("workflow", "modules", db_method + ".Rmd"),
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
		),
		samples = os.path.join(
			macs2_path, "{target}", "{target}_qc_samples.tsv"
		),
		scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
		setup = rules.create_setup_chunk.output,
		site_yaml = rules.create_site_yaml.output,
		yml = expand(
			os.path.join("config", "{file}.yml"), file = ['config', 'params']
		),
		rnaseq_mod = os.path.join(
			"workflow", "modules", "rnaseq_differential.Rmd"
		),
		windows = rules.count_windows.output
	output:
		html = "docs/{target}_{ref}_{treat}_differential_signal.html",
		fig_path = directory(
			os.path.join(
				"docs", "{target}_{ref}_{treat}_differential_signal_files",
				"figure-html"
			)
		),
		renv = temp(
			os.path.join(
				"output", "envs",
				"{target}_{ref}_{treat}-differential_signal.RData"
			)
		),
		outs = expand(
			os.path.join(
				diff_path, "{{target}}", "{{target}}_{{ref}}_{{treat}}-{f}"
			),
			f = [
				'differential_signal.rds', 'down.bed', 'up.bed','differential_signal.csv.gz', 'DE_genes.csv', 'enrichment.csv',
				'rnaseq_enrichment.csv'
			]
		)
	retries: 1
	conda: "../envs/rmarkdown.yml"
	threads:
		lambda wildcards: min(
			len(
				df[
					(df['target'] == wildcards.target) &
					(
						(df['treat'] == wildcards.ref) |
						(df['treat'] == wildcards.treat)
					)
				]
			),
			max_threads
		)
	resources:
		mem_mb = 65536,
		runtime = "3h"
	log: log_path + "/differential_signal/{target}_{ref}_{treat}_differential_signal.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""
