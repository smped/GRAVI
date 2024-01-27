rule create_annotations:
	input:
		bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = samples),
		blacklist = blacklist,
		config = ancient(os.path.join("config", "config.yml")),
		extrachips = rules.update_extrachips.output,
		gtf = gtf,
		r = os.path.join("workflow", "scripts", "create_annotations.R"),
		yaml = os.path.join("config", "params.yml")
	output:
		rds = expand(
		  os.path.join(annotation_path, "{file}.rds"),
		  file = [
			'gene_regions', 'gtf_gene', 'gtf_transcript', 'gtf_exon', 'seqinfo',
			'tss', 'trans_models'
			]
		),
		chrom_sizes = chrom_sizes
	params:
		annot_path = annotation_path
	conda: "../envs/rmarkdown.yml"
	threads: 1
	resources:
		mem_mb = 16384	
	log: log_path + "/scripts/create_annotations.log"
	shell:
		"""
		Rscript --vanilla {input.r} {params.annot_path} &>> {log}
		"""


rule compile_annotations_html:
    input:
        blacklist = blacklist,
        extrachips = rules.update_extrachips.output,
        here = here_file,
        rmd = "workflow/modules/annotation_description.Rmd",
        rds = rules.create_annotations.output.rds,
        scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
        setup = rules.create_setup_chunk.output,
        site_yaml = rules.create_site_yaml.output,
        yaml = expand(
            os.path.join("config", "{file}.yml"),
            file = ['config', 'colours', 'rmarkdown']
        )
    output:
        rmd = "analysis/annotation_description.Rmd",
        rds = os.path.join(annotation_path, "colours.rds"),
        html = "docs/annotation_description.html",
        fig_path = directory(
            os.path.join("docs", "annotation_description_files", "figure-html")
        )
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: log_path + "/rmarkdown/compile_annotations_html.log"
    resources:
        mem_mb = 4096,
        disk_mb = 4000,
    shell:
        """
        cp {input.rmd} {output.rmd}
        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}
        """
