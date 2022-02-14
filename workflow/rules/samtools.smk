rule index_bam:
	input:
		bam = os.path.join(bam_path, "{target}", "{sample}.bam")
	output:
		bai = os.path.join(bam_path, "{target}", "{sample}.bam.bai")
	conda: "../envs/samtools.yml"
	threads: 8
	shell:
		"""
		samtools index \
			-@ {threads} \
			{input.bam} \
			{output.bai}
		"""
