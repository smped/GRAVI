rule run_motif_analysis:
	input:
		exclude_ranges = os.path.join(annotation_path, "exclude_ranges.rds"),
		gene_regions = os.path.join(annotation_path, "gene_regions.rds"),
		motifs = os.path.join(annotation_path, "motif_list.rds"),
		packages = os.path.join(check_path, "r-packages.chk"),
		params = os.path.join("config", "params.yml"),
		peaks = os.path.join(
			peak_path, "{i}", "{i}_consensus_peaks.bed.gz"
		),
		script = os.path.join("workflow", "scripts", "motif_enrichment.R"),
	output:
		enrich = os.path.join(
			peak_path, "{i}", "{i}_motif_enrichment.tsv.gz"
		),
		pos = os.path.join(
			peak_path, "{i}", "{i}_motif_position.tsv.gz"
		),
	threads: 8
	retries: 1
	resources:
		runtime = "2h",
		mem_mb = 64000,
	log: os.path.join(log_path, "motif_analysis", "{i}_motif_enrichment.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/motif_enrichment.R"
