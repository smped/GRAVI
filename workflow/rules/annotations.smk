rule download_gtf:
	output: gtf
	params:
		url = urllib.parse.urlunparse(
			(
			'ftp', 'ftp.ebi.ac.uk',
			'pub/databases/gencode/Gencode_human/release_' +
			config['genome']['gencode'] + "/" + config['genome']['build'] +
			"_mapping/" + os.path.basename(gtf),
			'', '', ''
			)
		)
	log: "workflow/logs/downloads/download_gtf.log"
	shell:
		"""
		curl \
			-o {output} \
			{params.url} 2> {log}
		"""

rule download_blacklist:
	output: blacklist
	params:
		url = urllib.parse.urlunparse(
			(
				'https', 'raw.githubusercontent.com',
				'Boyle-Lab/Blacklist/master/lists/' +
				ucsc_build + "-blacklist.v2.bed.gz",
				'', '', ''
			)
		),
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	log: "workflow/logs/downloads/download_blacklist.log"
	shell:
		"""
		curl {params.url} --output {output} 2> {log}

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
			git add {output}
		fi
		"""

rule download_cytobands:
	output: cytobands
	params:
		url = urllib.parse.urlunparse(
			(
				'http', 'hgdownload.cse.ucsc.edu',
				'/goldenPath/' + ucsc_build + '/database/' +
				os.path.basename(cytobands),
				'', '', ''
			)
		),
		git = git_add
	log: "workflow/logs/downloads/download_cytobands.log"
	shell:
		"""
		curl \
			-o {output} \
			{params.url} 2> {log}
		if [[ {params.git} == "True" ]]; then
			git add {output}
		fi
		"""

rule write_chrom_sizes:
	input: cytobands
	output: chrom_sizes
	log: "workflow/logs/downloads/write_chrom_sizes.log"
	conda: "../envs/rmarkdown.yml"
	params:
		r = "workflow/scripts/write_chrom_sizes.R",
		git = git_add
	threads: 1
	shell:
		"""
		Rscript --vanilla \
			{params.r} \
			{input} \
			{output} &> {log}
		if [[ {params.git} == "True" ]]; then
			git add {output}
		fi
		"""


rule setup_annotations:
	input:
		blacklist = blacklist,
		bigwig = expand(
			os.path.join(bw_path, "{path}_merged_treat_pileup.bw"),
			path = merged_pre
		),
		chrom_sizes = chrom_sizes,
		config = "config/config.yml",
		gtf = gtf,
		## Figure how to make this optional later
		h3k27ac = expand(
			os.path.join(
				macs2_path, "H3K27ac", "{treat}_merged_peaks.narrowPeak"
			),
			treat = list(set(df[df.target == "H3K27ac"]['treat']))
		),
		pkgs = rules.install_packages.output,
		rnaseq = config['external']['rnaseq'],
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output,
		module = "workflow/modules/features_from_h3k27ac.Rmd",
		rmd = os.path.join(rmd_path, "annotation_setup.Rmd")
	output:
		annotations = ALL_ANNOTATION,
		html = "docs/annotation_setup.html",
		fig_path = directory(
			os.path.join("docs", "annotation_setup_files", "figure-html")
		),
		coverage = os.path.join(bw_path, "max_coverage.tsv")
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	threads: 8
	log: "workflow/logs/rmarkdown/annotation_setup.log"
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
			git add {output}
		fi
		"""

rule identify_super_enhancers:
	input:
		bam = expand(
			os.path.join(bam_path, "H3K27ac", "{sample}.bam"),
			sample = df[df.target == 'H3K27ac']['sample']
		),
		bai = expand(
			os.path.join(bam_path, "H3K27ac", "{sample}.bam.bai"),
			sample = df[df.target == 'H3K27ac']['sample']
		),
		bw = expand(
			os.path.join(
				bw_path, "H3K27ac", "{treat}_merged_treat_pileup.bw"
			),
			treat = list(set(df[df.target == "H3K27ac"]['treat']))
		),
		peaks = os.path.join("output", "H3K27ac", "consensus_peaks.bed"),
		annotations = ALL_ANNOTATION,
		blacklist = blacklist,
		coverage = rules.setup_annotations.output.coverage,
		cytobands = cytobands,
		pkgs = rules.install_packages.output,
		rnaseq = config['external']['rnaseq'],
		rose = "workflow/scripts/ROSE_callSuper.R",
		yaml = rules.create_site_yaml.output,
		setup = rules.create_setup_chunk.output,
		rmd = os.path.join(rmd_path, "rose.Rmd")
	output:
		html = "docs/rose.html",
		fig_path = directory("docs/rose_files/figure-html"),
		features = multiext(os.path.join(annotation_path, "SE"), ".rds", ".bed")
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	threads: len(df[df['target'] == "H3K27ac"])
	log: "workflow/logs/rmarkdown/rose.log"
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
			git add {output}
		fi
		"""
