rule create_annotations:
    input:
        bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = samples),
        chk = ALL_CHECKS,
        module = os.path.join(
            "workflow", "modules", "annotation_description.Rmd"
        ),
        script = os.path.join("workflow", "scripts", "create_annotations.R"),
        yaml = os.path.join("config", "params.yml"),
    output:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        exons = os.path.join(annotation_path, "gtf_exon.rds"),
        features = os.path.join(annotation_path, "features.rds"),
        genes = os.path.join(annotation_path, "gtf_gene.rds"),
        motifs = os.path.join(annotation_path, "motif_list.rds"),
        motif_uri = os.path.join(annotation_path, "motif_uri.rds"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        rmd = os.path.join(rmd_path, "annotation_description.Rmd"),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
        transcripts = os.path.join(annotation_path, "gtf_transcript.rds"),
        trans_models = os.path.join(annotation_path, "trans_models.rds"),
        tss = os.path.join(annotation_path, "tss.rds"),
        chrom_sizes = chrom_sizes
    conda: "../envs/rmarkdown.yml"
    params:
        annotation_path = annotation_path,
        colours = os.path.join(annotation_path, "colours.rds"),
        grey_path = grey_path,
    threads: 2
    resources:
        mem_mb = 16384,
        run_time = "30m"
    log: os.path.join(log_path, "scripts", "create_annotations.log")
    script:
        "../scripts/create_annotations.R"

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