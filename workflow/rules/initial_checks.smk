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
        
rule check_here_file:
    output: os.path.join("output", "checks", "here.chk")
    threads: 1
    resources:
        runtime = "1m",
        mem_mb = 1024,
    log: os.path.join(log_path, "checks", "check_here_file")
    shell:
        """
        f1=$(find ./ -type f -name '*Rproj')
        f2=$(find ./ -type f -name '*here')
        if [[ -z "$f1" && -z "$f2" ]]; then
          ## Create the file .here if neither .here nor Rproj exist
          echo "No viable here file detected. Creating .here"
          touch .here
          ## Check for success
          if [ -f "./.here" ]; then
            touch {output}
          fi
        else 
          echo "Found viable here file: $f1$f2"
          touch {output}
        fi
        """    