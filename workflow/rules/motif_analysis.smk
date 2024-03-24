## Is there a better way of handling 'i'?
## Maybe add 'shared' to the list of valid possible wildcards?
rule run_motif_analysis:
	input:
		exclude_ranges = os.path.join(annotation_path, "exclude_ranges.rds"),
		gene_regions = os.path.join(annotation_path, "gene_regions.rds"),
		motifs = os.path.join(annotation_path, "motif_list.rds"),
		packages = os.path.join(check_path, "r-packages.chk"),
		peaks = os.path.join(
			peak_path, "{target}", "{target}_consensus_peaks.rds"
		),
		script = os.path.join("workflow", "scripts", "motif_analysis.R"),
	output:
		enrich = os.path.join(
			peak_path, "{target}", "{target}_motif_enrichment.tsv.gz"
		),
		pos = os.path.join(
			peak_path, "{target}", "{target}_motif_position.tsv.gz"
		),
	params:
		abs =  lambda wildcards: motif_param[wildcards.target]['abs'],
		binwidth =  lambda wildcards: motif_param[wildcards.target]['binwidth'],
		ignore_below =  lambda wildcards: motif_param[wildcards.target]['ignore_below'],
		iterations =  lambda wildcards: motif_param[wildcards.target]['iterations'],
		model = lambda wildcards: motif_param[wildcards.target]['model'],
		peak_width = lambda wildcards: motif_param[wildcards.target]['peak_width'],
	threads: lambda wildcards, attempt: attempt * 8
	retries: 2
	resources:
		disk_mb = 10000,
		mem_mb = lambda wildcards, attempt: attempt * 64000,
		runtime = lambda wildcards, attempt: attempt * 120,
	log: os.path.join(log_path, "motif_analysis", "{target}_motif_analysis.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/motif_analysis.R"
