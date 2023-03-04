rule make_greylist:
	input: 
		bam = os.path.join(bam_path, "Input", "{ip_sample}.bam"),
		bai = os.path.join(bam_path, "Input", "{ip_sample}.bam.bai"),
		r = os.path.join("workflow", "scripts", "make_greylist.R"),
		seqinfo = os.path.join(annotation_path, "seqinfo.rds")
	output:
		bed = os.path.join(annotation_path, "{ip_sample}_greylist.bed")
	params:
		git = git_add,
	conda: "../envs/rmarkdown.yml"
	log: log_path + "/scripts/{ip_sample}_make_greylist.log"
	threads: 1
	retries: git_tries
	shell:
		"""
		Rscript --vanilla \
			{input.r} \
			{input.bam} \
			{input.seqinfo} \
			{output.bed} &>> {log}
		
		if [[ {params.git} == "True" ]]; then
			git add {output.bed}
		fi
		"""