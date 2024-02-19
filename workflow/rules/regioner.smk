rule regioner_localz_regions:
	input:
		checks = ALL_CHECKS,
		regions = os.path.join(annotation_path, "gene_regions.rds"),
		features = os.path.join(annotation_path, "features.rds"),
		peaks = os.path.join(
			macs2_path, "{path}_consensus_peaks.bed.gz"
		),
		params = os.path.join("config", "params.yml"),
	output:
		rds = os.path.join(macs2_path, "{path}_regions_localz.rds")
	params:
		ntimes = 5000, # Maybe source from elsewhere later?
	threads: 8
	resources:
		mem_mb = 32768,
		run_time = "60m",
	log: os.path.join(log_path, "regioner", "{path}_regions_localz.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/regioner_localz_regions.R"