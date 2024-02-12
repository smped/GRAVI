rule make_greylist:
	input: 
		bam = os.path.join(bam_path, "{ip_sample}.bam"),
		bai = os.path.join(bam_path, "{ip_sample}.bam.bai"),
		chk = ALL_CHECKS,
		r = os.path.join("workflow", "scripts", "make_greylist.R"),
		seqinfo = os.path.join(annotation_path, "seqinfo.rds")
	output:
		bed = os.path.join(annotation_path, "{ip_sample}_greylist.bed")
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/scripts/{ip_sample}_make_greylist.log"
	threads: 1
	resources:
		mem_mb = 16384
	script:
		"""
		{input.r}
		"""