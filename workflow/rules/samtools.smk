rule index_bam:
	input:
		bam = os.path.join(bam_path, "{bam}.bam")
	output:
		bai = os.path.join(bam_path, "{bam}.bam.bai")
	conda: "../envs/samtools.yml"
	threads: 8
	shell:
		"""
		samtools index -@ {threads} {input.bam} {output.bai}
		"""
