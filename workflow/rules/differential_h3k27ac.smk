
rule create_differential_h3k27ac_rmd:
	input:
		chk = rules.check_r_packages.output,
		db_mod = os.path.join(
			"workflow", "modules", "differential_h3k27ac.Rmd"
		),
		r = "workflow/scripts/create_differential_rmd.R"
	output: 
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_h3k27ac.Rmd"
		)
	params:
		git = git_add,
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
		type = "Signal"
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/create_rmd/{target}_{ref}_{treat}_differential_h3k27ac.log"
	threads: 1
	resources:
		runtime = "1m"
		mem_mb = 512
	retries: git_tries # Needed to subvert any issues with the git lock file
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

		if [[ {params.git} == "True" ]]; then
			git add {output.rmd}
		fi
		"""

rule count_h3k27ac_windows:
	input:
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{{target}}", "{sample}.{suffix}"),
			sample = df['sample'][df['target'] == wildcards.target],
			suffix = ['bam', 'bam.bai']
		),
		chk = rules.check_r_packages.output,
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
				"differential_h3k27ac", "{{target}}", "{{target}}_{file}.rds"
			),
			file = ['window_counts', 'filtered_counts']
		)
	params:
		type = "h3k27ac"
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/diferential_h3k27ac/{target}_count_h3k27ac_windows.log"
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


rule compile_differential_h3k27ac_html:
	input:
		annotations = ALL_RDS,
		chk = rules.check_r_packages.output,
		merged_macs2 = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", "{pre}_merged_callpeak.log"
			),
			pre = [wildcards.ref, wildcards.treat]
		),
		merged_bw = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", "{pre}_merged_treat_pileup.bw"
			),
			pre = [wildcards.ref, wildcards.treat]
		),
		peaks = os.path.join(macs2_path, "{target}", "consensus_peaks.bed"),
		here = here_file,
		indiv_bw = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", "{sample}_treat_pileup.bw"
			),
			sample = df['sample'][
				(df['target'] == wildcards.target) &
				(
					(df['treat'] == wildcards.ref) |
					(df['treat'] == wildcards.treat)
				)
			]
		),
		module = os.path.join("workflow", "modules", db_method + ".Rmd"),
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_h3k27ac.Rmd"
		),
		scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
		setup = rules.create_setup_chunk.output,
		site_yaml = rules.create_site_yaml.output,
		yml = expand(
			os.path.join("config", "{file}.yml"),
			file = ['config', 'params']
		),
		rnaseq_mod = os.path.join(
			"workflow", "modules", "rnaseq_differential.Rmd"
		),
		windows = rules.count_h3k27ac_windows.output
	output:
		html = "docs/{target}_{ref}_{treat}_differential_h3k27ac.html",
		fig_path = directory(
			os.path.join(
				"docs", "{target}_{ref}_{treat}_differential_h3k27ac_files",
				"figure-html"
			)
		),
		renv = temp(
			os.path.join(
				"output", "envs",
				"{target}_{ref}_{treat}-differential_h3k27ac.RData"
			)
		),
		outs = expand(
			os.path.join(
				diff_h3k27ac_path, 
				"{{target}}", "{{target}}_{{ref}}_{{treat}}-{f}"
			),
			f = [
				'differential_h3k27ac.rds', 'down.bed', 'up.bed','differential_h3k27ac.csv.gz', 'DE_genes.csv', 'enrichment.csv',
				'rnaseq_enrichment.csv'
			]
		)
	retries: 3
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = git_tries
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
		runtime = "1h",
		mem_mb = 16384		
	log: log_path + "/differential_h3k27ac/{target}_{ref}_{treat}_differential_h3k27ac.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
		"""
