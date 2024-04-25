rule check_r_packages:
    input: 
        here = os.path.join(check_path, "here.chk"),
        script = os.path.join("workflow", "scripts", "check_r_packages.R"),
        yml = os.path.join("workflow", "envs", "rmarkdown.yml"),
    output: os.path.join(check_path, "r-packages.chk")
    threads: 1
    resources:
        runtime = "30m",
        mem_mb = 2048,
    params:
        min_extrachips = "1.7.7"
    log: os.path.join(log_path, "initial_checks", "check_r_packages.log")
    conda: "../envs/rmarkdown.yml"
    localrule: True
    script:
        "../scripts/check_r_packages.R"
        
rule check_here_file:
    output: os.path.join(check_path, "here.chk")
    threads: 1
    localrule: True
    resources:
        runtime = "1m",
        mem_mb = 1024,
    log: os.path.join(log_path, "initial_checks", "check_here_file.log")
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

rule check_args:
  input:
    checks = rule.check_r_packages.output,
    colours = os.path.join("config", "colours.yml"),
    params = os.path.join("config", "params.yml"),
    script = os.path.join("workflow", "scripts", "check_all_args.R"),
  output: os.path.join(check_path, "args.chk")
  conda: "../envs/rmarkdown.yml"
  localrule: True
  resources:
    runtime = "10m",
    mem_mb = 1024,  
  script:
      "../scripts/check_all_args.R"
