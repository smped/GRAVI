rule localz_regions:
	input:
		checks = ALL_CHECKS,
		regions = os.path.join(annotation_path, "gene_regions.rds"),
		features = os.path.join(annotation_path, "features.rds"),
		peaks = os.path.join(peak_path, "{path}_consensus_peaks.bed.gz"),
		params = os.path.join("config", "params.yml"),
	output:
		rds = os.path.join(peak_path, "{path}_regions_localz.rds")
	threads: 8
	retries: 1
	resources:
		mem_mb = 32768,
		run_time = "60m",
	log: os.path.join(log_path, "regioner", "{path}_regions_localz.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/regioner_localz_regions.R"

rule localz_targets:
	input:
		checks = ALL_CHECKS,
		params = os.path.join("config", "params.yml"),
		peaks = expand(
			os.path.join(
				peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
			),
			target = targets
		),
		sq = os.path.join(annotation_path, "seqinfo.rds")
	output:
		rds = os.path.join(peak_path, "shared", "shared_targets_localz.rds")
	threads: 8
	retries: 1
	resources:
		mem_mb = 32768,
		run_time = "60m",
	log: os.path.join(log_path, "regioner", "shared", "shared_targets_localz.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/regioner_localz_targets.R"