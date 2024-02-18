rule check_r_packages:
    input: 
        script = os.path.join("workflow", "scripts", "check_r_packages.R"),
        yml = os.path.join("workflow", "envs", "rmarkdown.yml"),
    output: os.path.join("output", "checks", "r-packages.chk")
    threads: 1
    resources:
        runtime = "30m",
        mem_mb = 2048,
    log: os.path.join(log_path, "checks", "check_r_packages.log")
    conda: "../envs/rmarkdown.yml"
    localrule: True
    script:
        "../scripts/check_r_packages.R"
        
rule check_here_file:
    output: os.path.join("output", "checks", "here.chk")
    threads: 1
    localrule: True
    resources:
        runtime = "1m",
        mem_mb = 1024,
    log: os.path.join(log_path, "checks", "check_here_file.log")
    shell:
        """
        f1=$(find ./ -type f -name '*Rproj')
        f2=$(find ./ -type f -name '*here')
        if [[ -z "$f1" && -z "$f2" ]]; then
          ## Create the file .here if neither .here nor Rproj exist
          echo "No viable here file detected. Creating .here" >> {log}
          touch .here
          ## Check for success
          if [ -f "./.here" ]; then
            touch {output}
          fi
        else 
          echo "Found viable here file: $f1$f2" >> {log}
          touch {output}
        fi
        """    

rule check_external_files:
    input: 
        bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = samples),
        bai = expand(os.path.join(bam_path, "{bam}.bam.bai"), bam = samples),
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        script = "../scripts/check_external_files.R",
    output: os.path.join("output", "checks", "external-files.chk")
    threads: 1
    resources:
        runtime = "30m",
        mem_mb = 8192,
    log: os.path.join(log_path, "checks", "check_external_files.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/check_external_files.R"