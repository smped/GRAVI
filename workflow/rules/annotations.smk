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

rule setup_annotations:
	input:
		bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = indiv_pre),
		blacklist = blacklist,
		config = "config/config.yml",
		gtf = gtf,
		pkgs = rules.install_packages.output,
		rnaseq = config['external']['rnaseq'],
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output,
		rmd = os.path.join(rmd_path, "annotation_setup.Rmd")
	output:
		annotations = ALL_ANNOTATION,
		chrom_sizes = chrom_sizes,
		html = "docs/annotation_setup.html",
		fig_path = directory(
			os.path.join("docs", "annotation_setup_files", "figure-html")
		)
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

