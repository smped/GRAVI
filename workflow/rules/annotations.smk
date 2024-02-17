rule create_annotations:
    input:
        bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = samples),
        chk = ALL_CHECKS,
        yaml = os.path.join("config", "params.yml"),
    output:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        exons = os.path.join(annotation_path, "gtf_exon.rds"),
        features = os.path.join(annotation_path, "features.rds"),
        genes = os.path.join(annotation_path, "gtf_gene.rds"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
        transcripts = os.path.join(annotation_path, "gtf_transcript.rds"),
        trans_models = os.path.join(annotation_path, "trans_models.rds"),
        tss = os.path.join(annotation_path, "tss.rds"),
        chrom_sizes = chrom_sizes
    conda: "../envs/rmarkdown.yml"
    threads: 2
    resources:
        mem_mb = 16384,
        run_time = "30m"
    log: os.path.join(log_path, "scripts", "create_annotations.log")
    script:
        "../scripts/create_annotations.R"

rule compile_annotations_html:
    input:
        checks = ALL_CHECKS,
        greylist = expand(
            os.path.join(annotation_path, "{f}_greylist.bed.gz"),
            f = set(df['input'])
        ),
        rmd = "workflow/modules/annotation_description.Rmd",
        rds = rules.create_annotations.output,
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
    log: os.path.join(log_path, "rmarkdown", "compile_annotations_html.log")
    resources:
        mem_mb = 4096,
        disk_mb = 4000,
        run_time = "10m",
    shell:
        """
        cp {input.rmd} {output.rmd}
        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}
        """

