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
		type = "Binding"
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/create_rmd/{target}_{ref}_{treat}_differential_binding.log"
	threads: 1
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

rule count_tf_windows:
	input:
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
			len(df[(df['target'] == wildcards.target)]),
			max_threads
		)
	shell:
		"""
		Rscript --vanilla \
		  {input.r} \
		  {wildcards.targets} \
		  {threads} 
		  {params.type} &>> {log}
		"""

rule compile_differential_binding_html:
	input:
		annotations = ALL_RDS,
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{{target}}", "{sample}.{suffix}"),
			sample = df['sample'][
				(df['target'] == wildcards.target) &
				(
					(df['treat'] == wildcards.ref) |
					(df['treat'] == wildcards.treat)
				)
			],
			suffix = ['bam', 'bam.bai']
		),
		extrachips = rules.update_extrachips.output,
		greylist = lambda wildcards: expand(
			os.path.join(annotation_path, "{ip_sample}_greylist.bed"),
			ip_sample = set(
				df['input'][
					(df['target'] == wildcards.target) &
					(
						(df['treat'] == wildcards.ref) |
						(df['treat'] == wildcards.treat)
					)
				]
			),
		),
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
			rmd_path, "{target}_{ref}_{treat}_differential_binding.Rmd"
		),
		samples = os.path.join(macs2_path, "{target}", "qc_samples.tsv"),
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
		windows = os.path.join(
			diff_tf_path, "{target}", "{target}_filtered_counts.rds"
		)
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
	log: log_path + "/differential_binding/{target}_{ref}_{treat}_differential_binding.log"
	shell:
		"""
		R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}

		if [[ {params.git} == "True" ]]; then
			TRIES={params.tries}
			while [[ -f .git/index.lock ]]
			do
				if [[ "$TRIES" == 0 ]]; then
					echo "ERROR: Timeout while waiting for removal of git index.lock" &>> {log}
					exit 1
				fi
				sleep {params.interval}
				((TRIES--))
			done
			git add {output.html} {output.outs}
			git add {output.fig_path}
		fi
		"""
