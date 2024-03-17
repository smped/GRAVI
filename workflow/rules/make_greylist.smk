rule make_greylist:
    input: 
        bam = os.path.join(bam_path, "{ip_sample}.bam"),
        bai = os.path.join(bam_path, "{ip_sample}.bam.bai"),
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        sq = os.path.join(annotation_path, "seqinfo.rds")
    output:
        bed = os.path.join(grey_path, "{ip_sample}_greylist.bed.gz")
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "make_greylist", "{ip_sample}_make_greylist.log")
    threads: 2
    resources:
        mem_mb = 16384,
        run_time = "30m"
    script:
        "../scripts/make_greylist.R"
