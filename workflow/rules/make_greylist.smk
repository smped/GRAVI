rule make_greylist:
    input: 
        bam = os.path.join(bam_path, "{ip_sample}.bam"),
        bai = os.path.join(bam_path, "{ip_sample}.bam.bai"),
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "make_greylist.R"),
        sq = rules.create_genome_annotations.output.seqinfo, 
    output:
        bed = os.path.join(grey_path, "{ip_sample}_greylist.bed.gz")
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "greylist", "{ip_sample}_make_greylist.log")
    threads: 2
    retries: 1    
    resources:
        mem_mb = 16384,
        runtime = "30m"
    script:
        "../scripts/make_greylist.R"

rule combine_greylists:
    input:
        gl = expand(
            os.path.join(grey_path, "{f}_greylist.bed.gz"),
            f = set(df['input'])
        ),
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "combine_greylists.R"),
        sq = rules.create_genome_annotations.output.seqinfo, 
    output:
        rds = os.path.join(grey_path, "greylists.rds")
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "greylist", "combine_greylists.log")
    threads: 2
    retries: 1    
    resources:
        mem_mb = 16000,
        runtime = "20m"
    script:
        "../scripts/combine_greylists.R"    