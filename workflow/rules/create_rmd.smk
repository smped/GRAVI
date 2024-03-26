rule create_site_yaml:
    input:
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "create_site_yaml.R"),
        yml = "config/rmarkdown.yml",
    output:
        yml = os.path.join(rmd_path, "_site.yml")
    log: os.path.join(log_path, "scripts", "create_site_yaml.log")
    params:
        targets = targets,
        diff_sig = diff_sig_param,
        pairs = pairs,
    threads: 1
    localrule: True
    resources:
        mem_mb = 1024,
        runtime = "5m",
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/create_site_yaml.R"

rule create_setup_chunk:
    input:
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "create_setup_chunk.R"),
        yml = "config/rmarkdown.yml",
    output:
        rmd = "analysis/setup_chunk.Rmd"
    log: os.path.join(log_path, "scripts", "create_setup_chunk.log")
    threads: 1
    localrule: True
    resources:
        mem_mb = 1024,
        runtime = "5m",
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/create_setup_chunk.R"

rule create_index_rmd:
    input:
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        rmd = os.path.join("workflow", "modules", "index.Rmd"),
    output:
        os.path.join(rmd_path, "index.Rmd")
    threads: 1
    localrule: True
    resources:
        mem_mb = 512,
        runtime = "2m",
    shell:
        """
        cat {input.rmd} > {output}
        """

rule create_macs2_summary_rmd:
    input:
        here = rules.check_here_file.output,
        module = "workflow/modules/macs2_summary.Rmd",
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "create_macs2_summary.R"),
    output:
        rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd")
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "create_rmd", "create_{target}_macs2_summary.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_macs2_summary.R"

rule create_differential_signal_rmd:
    input:
        chk = ALL_CHECKS,
        counts = os.path.join(diff_path, "{target}", "{target}_counts.rds"),
        module = os.path.join("workflow", "modules", "differential_signal.Rmd"),
        r = os.path.join("workflow", "scripts", "create_differential_rmd.R")
    output:
        rmd = os.path.join(
            rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
        )
    params:
        diff_sig_params = lambda wildcards: diff_sig_param[wildcards.target]
    conda: "../envs/rmarkdown.yml"
    localrule: True
    log: os.path.join(log_path, "create_rmd", "{target}_{ref}_{treat}_differential_signal.log")
    threads: 1
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_differential_rmd.R"

