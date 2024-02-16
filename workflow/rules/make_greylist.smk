rule make_greylist:
	input: 
		bam = os.path.join(bam_path, "{ip_sample}.bam"),
		bai = os.path.join(bam_path, "{ip_sample}.bam.bai"),
		chk = ALL_CHECKS,
		sq = os.path.join(annotation_path, "seqinfo.rds")
	output:
		bed = os.path.join(annotation_path, "{ip_sample}_greylist.bed.gz")
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/scripts/{ip_sample}_make_greylist.log"
	threads: 1
	resources:
		mem_mb = 16384
	script:
		"../scripts/make_greylist.R"