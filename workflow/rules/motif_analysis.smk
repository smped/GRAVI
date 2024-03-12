rule make_exclude_ranges:
	input:
		here = os.path.join(check_path, "here.chk"),
		packages = os.path.join(check_path, "r-packages.chk"),
		script = os.path.join(
			"workflow", "scripts", "make_exclude_ranges.R"
		),
		seqinfo = os.path.join(annotation_path, "seqinfo.rds")
	output:
		rds = os.path.join(annotation_path, "exclude_ranges.rds")
	threads: 1
	localrule: True
	resources:
		runtime = "10m",
		mem_mb = 4096,
	log: os.path.join(log_path, "scripts", "make_exclude_ranges.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/make_exclude_ranges.R"

rule run_motif_analysis:
	input:
		exclude_ranges = os.path.join(annotation_path, "exclude_ranges.rds"),
		gene_regions = os.path.join(annotation_path, "gene_regions.rds"),
		motifs = os.path.join(annotation_path, "motif_list.rds"),
		packages = os.path.join(check_path, "r-packages.chk"),
		params = os.path.join("config", "params.yml"),
		peaks = os.path.join(
			macs2_path, "{i}", "{i}_consensus_peaks.bed.gz"
		),
		script = os.path.join("workflow", "scripts", "motif_enrichment.R"),
	output:
		enrich = os.path.join(
			macs2_path, "{i}", "{i}_motif_enrichment.tsv.gz"
		),
		matches = os.path.join(macs2_path, "{i}", "{i}_matches.rds"),
		pos = os.path.join(
			macs2_path, "{i}", "{i}_motif_position.tsv.gz"
		),
	threads: 8
	resources:
		runtime = "2h",
		mem_mb = 32768,
	log: os.path.join(log_path, "motif_analysis", "{i}_motif_enrichment.log")
	conda: "../envs/rmarkdown.yml"
	script:
		"../scripts/motif_enrichment.R"
