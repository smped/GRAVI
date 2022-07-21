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
	log: log_path + "/downloads/download_gtf.log"
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
	log: log_path + "/downloads/download_blacklist.log"
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

rule create_annotations:
	input:
		bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = indiv_pre),
		config = ancient(os.path.join("config", "config.yml")),
		gtf = gtf,
		r = os.path.join("workflow", "scripts", "create_annotations.R"),
		pkgs = rules.install_packages.output,
		yaml = os.path.join("config", "params.yml")
	output:
		rds = expand(
		  os.path.join(annotation_path, "{file}.rds"),
		  file = ['all_gr', 'gene_regions', 'tss', 'seqinfo', 'trans_models']
		),
		chrom_sizes = chrom_sizes
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = git_tries
	conda: "../envs/rmarkdown.yml"
	threads: max_threads
	log: log_path + "/scripts/create_annotations.log"
	shell:
		"""
		Rscript --vanilla {input.r} {input.gtf} {threads} &>> {log}

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

