rule check_r_packages:
    input: os.path.join("workflow", "scripts", "check_r_packages.R")
    output: os.path.join("output", "checks", "r-packages.chk")
    conda: "../envs/rmarkdown.yml"
    threads: 1
    resources:
        runtime = "30m",
        mem_mb = 4096,
    log: os.path.join(log_path, "checks", "check_r_packages")
    shell:
        """
        Rscript --vanilla {input} {output} &>> {log}
        """
        