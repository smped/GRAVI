rule make_greylist:
	input: 
		bam = os.path.join(bam_path, "{ip_sample}.bam"),
		bim = os.path.join(bam_path, "{ip_sample}.bam.bai"),
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

rule create_differential_binding_rmd:
	input:
		db_mod = os.path.join(
			"workflow", "modules", "differential_binding.Rmd"
		),
		r = "workflow/scripts/create_differential_rmd.R"
	output: 
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_binding.Rmd"
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
		type = "Binding"
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/create_rmd/{target}_{ref}_{treat}_differential_binding.log"
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
			{params.type} \
			{output.rmd} &>> {log}

		## Add the remainder of the module as literal text
		cat {input.db_mod} >> {output.rmd}
		"""

rule count_tf_windows:
	input:
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{sample}.{suffix}"),
			sample = df['sample'][
				(df['target'] == wildcards.target) &
				(
					(df['treat'] == wildcards.ref) |
					(df['treat'] == wildcards.treat)
				)
			],
			suffix = ['bam', 'bam.bai']
		),
		blacklist = blacklist,
		extrachips = rules.update_extrachips.output,
		greylist = lambda wildcards: expand(
			os.path.join(annotation_path, "{ip_sample}_greylist.bed"),
			ip_sample = set(df['input'][df['target'] == wildcards.target])
		),
		peaks = os.path.join(macs2_path, "{target}", "consensus_peaks.bed"),
		samples = os.path.join(macs2_path, "{target}", "qc_samples.tsv"),
		r = os.path.join("workflow", "scripts", "make_filtered_counts.R")
	output:
		rds = expand(
			os.path.join(
				"differential_binding", "{{target}}", "{{target}}_{file}.rds"
			),
			file = ['window_counts', 'filtered_counts']
		)
	params:
		type = "binding"
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/diferential_binding/{target}_count_tf_windows.log"
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
		  {threads} 
		  {params.type} &>> {log}
		"""

rule compile_differential_binding_html:
	input:
		annotations = ALL_RDS,
		extrachips = rules.update_extrachips.output,
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
		here = here_file,
		module = os.path.join("workflow", "modules", db_method + ".Rmd"),
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_binding.Rmd"
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
		windows = rules.count_tf_windows.output
	output:
		html = "docs/{target}_{ref}_{treat}_differential_binding.html",
		fig_path = directory(
			os.path.join(
				"docs", "{target}_{ref}_{treat}_differential_binding_files",
				"figure-html"
			)
		),
		renv = temp(
			os.path.join(
				"output", "envs",
				"{target}_{ref}_{treat}-differential_binding.RData"
			)
		),
		outs = expand(
			os.path.join(
				diff_tf_path, "{{target}}", "{{target}}_{{ref}}_{{treat}}-{f}"
			),
			f = [
				'differential_binding.rds', 'down.bed', 'up.bed','differential_binding.csv.gz', 'DE_genes.csv', 'enrichment.csv',
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
	log: log_path + "/differential_binding/{target}_{ref}_{treat}_differential_binding.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""
