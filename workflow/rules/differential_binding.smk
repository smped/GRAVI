rule create_differential_binding_rmd:
	input:
		db_mod = os.path.join(
			"workflow", "modules", "differential_binding.Rmd"
		),
		r = "workflow/scripts/create_differential_binding.R"
	output: 
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_binding.Rmd"
		)
	params:
		git = git_add,
		threads = lambda wildcards: len(df[
			(df['target'] == wildcards.target) &
			((df['treat'] == wildcards.ref) | (df['treat'] == wildcards.treat))
		])
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
			{output.rmd} &>> {log}

		## Add the remainder of the module as literal text
		cat {input.db_mod} >> {output.rmd}

		if [[ {params.git} == "True" ]]; then
			git add {output.rmd}
		fi
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
		peaks = expand(
			os.path.join(macs2_path, "{target}", "consensus_peaks.bed"),
			target = targets
		),
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
		pkgs = rules.install_packages.output,
		rmd = os.path.join(
			rmd_path, "{target}_{ref}_{treat}_differential_binding.Rmd"
		),
		samples = os.path.join(macs2_path, "{target}", "qc_samples.tsv"),
		setup = rules.create_setup_chunk.output,
		site_yaml = rules.create_site_yaml.output,
		yml = expand(
			os.path.join("config", "{file}.yml"),
			file = ['config', 'params']
		),
		rnaseq_mod = os.path.join(
			"workflow", "modules", "rnaseq_differential_binding.Rmd"
		)
	output:
		html = "docs/{target}_{ref}_{treat}_differential_binding.html",
		fig_path = directory(
			os.path.join(
				"docs", "{target}_{ref}_{treat}_differential_binding_files",
				"figure-html"
			)
		),
		renv = expand(
			os.path.join(
				"output", "envs",
				"{{target}}_{{ref}}_{{treat}}-differential_binding.RData"
			)
		),
		outs = expand(
			os.path.join(
				diff_path, "{{target}}", "{{target}}_{{ref}}_{{treat}}-{file}"
			),
			file = [
				'differential_binding.rds', 'down.bed', 'up.bed','differential_binding.csv.gz', 'DE_genes.csv'
			]
		),
		win = os.path.join(
			diff_path, "{target}", "{target}_{ref}_{treat}-filtered_windows.rds"
		)
	retries: 3
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = git_tries
	conda: "../envs/rmarkdown.yml"
	threads:
		lambda wildcards: len(df[
			(df['target'] == wildcards.target) &
			((df['treat'] == wildcards.ref) | (df['treat'] == wildcards.treat))
			])
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
			git add {output.html}
			git add {output.fig_path}
			git add {output.outs}
		fi
		"""
