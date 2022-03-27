def get_input_bam_from_sample_and_target(wildcards):
	ind = (df['sample'] == wildcards.sample) & (df['target'] == wildcards.target)
	return expand(
		os.path.join(bam_path, "Input", "{file}.bam"),
		file = set(df[ind]['input'])
	)

def get_input_bai_from_sample_and_target(wildcards):
	ind = (df['sample'] == wildcards.sample) & (df['target'] == wildcards.target)
	return expand(
		os.path.join(bam_path, "Input", "{file}.bam.bai"),
		file = set(df[ind]['input'])
	)

def get_merged_bam_from_treat_and_target(wildcards):
	ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
	return expand(
		os.path.join(bam_path, wildcards.target, "{file}.bam"),
		file = set(df[ind]['sample'])
	)

def get_merged_bai_from_treat_and_target(wildcards):
	ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
	return expand(
		os.path.join(bam_path, wildcards.target, "{file}.bam.bai"),
		file = set(df[ind]['sample'])
	)

def get_input_bam_from_treat_and_target(wildcards):
	ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
	return expand(
		os.path.join(bam_path, "Input", "{file}.bam"),
		file = set(df[ind]['input'])
	)

def get_input_bai_from_treat_and_target(wildcards):
	ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
	return expand(
		os.path.join(bam_path, "Input", "{file}.bam.bai"),
		file = set(df[ind]['input'])
	)

rule macs2_individual:
	input:
		bam = os.path.join(bam_path, "{target}", "{sample}.bam"),
		bai = os.path.join(bam_path, "{target}", "{sample}.bam.bai"),
		control = get_input_bam_from_sample_and_target,
		control_bai = get_input_bai_from_sample_and_target
	output:
		narrow_peaks = os.path.join(
			macs2_path, "{target}", "{sample}_peaks.narrowPeak"
		),
		bedgraph = temp(
			os.path.join(
				macs2_path, "{target}", "{sample}_treat_pileup.bdg"
			)
		),
		log = os.path.join(
			macs2_path, "{target}", "{sample}_callpeak.log"
		),
		other = temp(
			expand(
				os.path.join(
					macs2_path, "{{target}}", "{{sample}}{suffix}"
				),
				suffix = [
					'_model.r', '_peaks.xls', '_summits.bed',
					'_control_lambda.bdg'
				]
			)
		)
	log: "workflow/logs/macs2_individual/{target}/{sample}.log"
	conda: "../envs/macs2.yml"
	params:
		outdir = os.path.join(macs2_path, "{target}"),
		prefix = "{sample}",
		gsize = config['peaks']['macs2']['gsize'],
		fdr = config['peaks']['macs2']['fdr'],
		keep_duplicates = config['peaks']['macs2']['keep_duplicates'],
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	threads: 1
	shell:
		"""
		echo -e "Running macs2 call peak on:\n{input.bam}" >> {log}
		echo -e "The specified control sample is:\n{input.control}" >> {log}
		macs2 callpeak \
			-t {input.bam}\
			-c {input.control} \
			-f BAM \
			-g {params.gsize} \
			--keep-dup {params.keep_duplicates} \
			-q {params.fdr} \
			-n {params.prefix} \
			--bdg --SPMR \
			--outdir {params.outdir} 2> {output.log}

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
			git add {output.narrow_peaks}
			git add {output.log}
		fi
		"""

rule macs2_qc:
	input:
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{{target}}", "{sample}.{suffix}"),
			sample = set(df[df.target == wildcards.target]['sample']),
			suffix = ['bam', 'bam.bai']
		),
		blacklist = blacklist,
		config = "config/config.yml",
		indiv_macs2 = lambda wildcards: expand(
			os.path.join(macs2_path, "{{target}}", "{sample}_{suffix}"),
			sample = set(df[df.target == wildcards.target]['sample']),
			suffix = ['callpeak.log', 'peaks.narrowPeak']
		),
		pkgs = rules.install_packages.output,
		samples = config['samples']['file'],
		seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
		r = "workflow/scripts/macs2_qc.R"
	output:
		cors = os.path.join("output", "{target}", "cross_correlations.tsv"),
		qc = os.path.join("output", "{target}", "qc_samples.tsv")
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	threads: lambda wildcards: len(df[df['target'] == wildcards.target])
	log: "workflow/logs/macs2_individual/{target}/{target}_macs2_qc.log"
	shell:
		"""
		## Run the QC script
		Rscript --vanilla \
			{input.r} \
			{wildcards.target} \
			{threads} &>> {log}

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

rule macs2_merged:
	input:
		bam = get_merged_bam_from_treat_and_target,
		bai = get_merged_bai_from_treat_and_target,
		control = get_input_bam_from_treat_and_target,
		control_bai = get_input_bai_from_treat_and_target,
		qc = os.path.join("output", "{target}", "qc_samples.tsv")
	output:
		narrow_peaks = os.path.join(
			macs2_path, "{target}", "{treat}_merged_peaks.narrowPeak"
		),
		bedgraph = temp(
			os.path.join(
				macs2_path, "{target}", "{treat}_merged_treat_pileup.bdg"
			)
		),
		log = os.path.join(
			macs2_path, "{target}", "{treat}_merged_callpeak.log"
		),
	log: "workflow/logs/macs2_merged/{target}/{treat}_merged.log"
	conda: "../envs/macs2.yml"
	params:
		bamdir = bam_path,
		outdir = os.path.join(macs2_path, "{target}"),
		prefix = "{treat}_merged",
		gsize = config['peaks']['macs2']['gsize'],
		fdr = config['peaks']['macs2']['fdr'],
		keep_duplicates = config['peaks']['macs2']['keep_duplicates'],
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	threads: 1
	shell:
		"""
		QC_PASS=$(egrep 'pass' {input.qc} | egrep {wildcards.treat} | cut -f1)
		SAMPLES=$(for f in $QC_PASS; do echo {params.bamdir}/{wildcards.target}/$f.bam; done)
		echo -e "The following passed qc:\n$SAMPLES\n" >> {log}
		INPUT_PASS=$(egrep 'pass$' {input.qc} | egrep 'DHT' | cut -f6 | uniq)
		INPUT=$(for f in $INPUT_PASS; do echo {params.bamdir}/Input/$f.bam; done)
		echo -e "Control:\n$INPUT" >> {log}
		macs2 callpeak \
			-t $SAMPLES\
			-c $INPUT \
			-f BAM \
			-g {params.gsize} \
			--keep-dup {params.keep_duplicates} \
			-q {params.fdr} \
			-n {params.prefix} \
			--bdg --SPMR \
			--outdir {params.outdir} 2> {output.log}

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
			git add {output.narrow_peaks}
			git add {output.log}
		fi
		"""

rule bedgraph_to_bigwig:
	input:
		bedgraph = os.path.join(
			macs2_path, "{target}", "{sample}_treat_pileup.bdg"
		),
		chrom_sizes = chrom_sizes
	output:
		bigwig = os.path.join(bw_path, "{target}", "{sample}_treat_pileup.bw")
	conda: "../envs/bedgraph_to_bigwig.yml"
	log: "workflow/logs/bedgraph_to_bigwig/{target}/{sample}.log"
	threads: 1
	shell:
		"""
		echo -e "\nConverting {input.bedgraph} to BigWig\n" >> {log}
		TEMPDIR=$(mktemp -d -t bdgXXXXXXXXXX)
		SORTED_BDG=$TEMPDIR/temp.bdg

		## Sort the file
		echo -e "\nSorting as $SORTED_BDG..." >> {log}
		sort -k1,1 -k2,2n {input.bedgraph} | \
			egrep '^chr[0-9XY]' > $SORTED_BDG

		## Convert the file
		echo -e "Done\nConverting..." >> {log}
		bedGraphToBigWig \
			$SORTED_BDG \
			{input.chrom_sizes} \
			{output.bigwig}
		echo -e "Done" >> {log}

		## Remove the temp sorted file
		rm -rf $TEMPDIR
		"""

rule get_coverage_summary:
	input: rules.bedgraph_to_bigwig.output.bigwig
	output: os.path.join(bw_path, "{target}", "{sample}_treat_pileup.summary")
	params:
		script = "workflow/scripts/get_bigwig_summary.R",
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	log: "workflow/logs/get_coverage_summary/{target}/{sample}.log"
	threads: 1
	shell:
		"""
		## Create the summary tsv
		Rscript --vanilla \
			{params.script} \
			{input} \
			{output} &>> {log}

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

rule build_macs2_summary:
	input:
		annotations = ALL_RDS,
		aln = lambda wildcards: expand(
			os.path.join(bam_path, "{{target}}", "{sample}.{suffix}"),
			sample = set(df[df.target == wildcards.target]['sample']),
			suffix = ['bam', 'bam.bai']
		),
		blacklist = blacklist,
		indiv_macs2 = lambda wildcards: expand(
			os.path.join(macs2_path, "{{target}}", "{sample}_{suffix}"),
			sample = set(df[df.target == wildcards.target]['sample']),
			suffix = ['callpeak.log', 'peaks.narrowPeak']
		),
		merged_macs2 = lambda wildcards: expand(
			os.path.join(
				macs2_path, "{{target}}", "{treat}_merged_{suffix}"
			),
			treat = set(df[df.target == wildcards.target]['treat']),
			suffix = ['callpeak.log', 'peaks.narrowPeak']
		),
		cors = os.path.join("output", "{target}", "cross_correlations.tsv"),
		qc = os.path.join("output", "{target}", "qc_samples.tsv"),
		config = "config/config.yml",
		pkgs = rules.install_packages.output,
		r = "workflow/scripts/create_macs2_summary.R",
		setup = rules.create_setup_chunk.output,
		yaml = rules.create_site_yaml.output,
		module = "workflow/modules/macs2_summary.Rmd"
	output:
		rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd"),
		html = "docs/{target}_macs2_summary.html",
		fig_path = directory(
			os.path.join(
				"docs", "{target}_macs2_summary_files",
				"figure-html"
			)
		),
		venn = "docs/assets/{target}/{target}_common_peaks.png",
		peaks = expand(
			"output/{{target}}/{file}",
			file = ['consensus_peaks.bed', 'oracle_peaks.rds']
		)
	params:
		git = git_add,
		interval = random.uniform(0, 1),
		tries = 10
	conda: "../envs/rmarkdown.yml"
	threads: lambda wildcards: len(df[df['target'] == wildcards.target])
	log: "workflow/logs/rmarkdown/build_{target}_macs2_summary.log"
	shell:
		"""
		## Create the generic markdown
		Rscript --vanilla \
			{input.r} \
			{wildcards.target} \
			{threads} \
			{output.rmd} &>> {log}

		R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}

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
