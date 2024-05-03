rule make_greylist:
    input: 
        bam = os.path.join(bam_path, "{ip_sample}.bam"),
        bai = os.path.join(bam_path, "{ip_sample}.bam.bai"),
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "make_greylist.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds")
    output:
        bed = os.path.join(grey_path, "{ip_sample}_greylist.bed.gz")
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "greylist", "{ip_sample}_make_greylist.log")
    threads: 2
    resources:
        mem_mb = 16384,
        run_time = "30m"
    script:
        "../scripts/make_greylist.R"

rule combine_greylists:
    input:
        gl = expand(
            os.path.join(grey_path, "{f}_greylist.bed.gz"),
            f = [set(df['input'])]
        ),
        script = os.path.join("workflow", "scripts", "combine_greylists.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds")
    output:
        rds = os.path.join(grey_path, "greylists.rds")
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "greylist", "combine_greylists.log")
    threads: 2
    resources:
        mem_mb = 16000,
        run_time = "20m"
    script:
        "../scripts/combine_greylists.R"    