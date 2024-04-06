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

rule create_signal_summary_rmd:
    input:
        annotations = ANNOTATION_RDS,
        bw = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}",
                "{{target}}_{treat}_merged_treat_pileup.bw"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        cors = os.path.join(
            macs2_path, "{target}", "{target}_cross_correlations.tsv"
        ),
        external = rules.check_external_files.output,
        here = rules.check_here_file.output,
        module = "workflow/modules/signal_summary.Rmd",
        packages = rules.check_r_packages.output,
        peak_files = expand(
             os.path.join(
                peak_path, "{{target}}", "{{target}}_{f}"
            ),
            f = [
                'motif_position.tsv.gz', 'motif_enrichment.tsv.gz',
                'regions_localz.rds', 'consensus_peaks.bed.gz'
                ]
        ),
        script = os.path.join("workflow", "scripts", "create_signal_summary.R"),
    output:
        rmd = os.path.join(rmd_path, "{target}_signal_summary.Rmd")
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "create_rmd", "create_{target}_signal_summary.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_signal_summary.R"


rule create_nfr_rmd:
    input:
        consensus_peaks = expand(
            os.path.join(peak_path, "{t}", "{t}_consensus_peaks.rds"),
            t = targets
        ),
        here = rules.check_here_file.output,
        files = expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_consensus_nfr.{suffix}"
            ),
            suffix = ['rds', 'bed.gz']
        ),
        localz = expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_nfr_{f}_localz.rds"
            ),
            f = ['regions', 'targets']
        ),
        module = os.path.join("workflow", "modules", "nfr.Rmd"),
        motif_results = expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_motif_{f}.tsv.gz"
            ),
            f = ['enrichment', 'position']
        ),
        nfr = lambda wildcards: expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_{treat}.nfr.bed.gz"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        packages = rules.check_r_packages.output,
        peaks = lambda wildcards: expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_{treat}_merged_peaks.bed"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        script = os.path.join("workflow", "scripts", "create_nfr_rmd.R"),
    output:
        rmd = os.path.join(rmd_path, "{target}_nfr.Rmd")
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "create_rmd", "create_{target}_nfr.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_nfr_rmd.R"

rule create_differential_signal_rmd:
    input:
        annotations = ANNOTATION_RDS,
        chk = ALL_CHECKS,
        counts = os.path.join(diff_path, "{target}", "{target}_counts.rds"),
        ihw = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-ihw.rds"
        ),
        results = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}-differential-signal.rds"
        ),
        module = os.path.join("workflow", "modules", "differential_signal.Rmd"),
        r = os.path.join("workflow", "scripts", "create_differential_rmd.R")
    output:
        rmd = os.path.join(
            rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
        )
    conda: "../envs/rmarkdown.yml"
    localrule: True
    log: os.path.join(log_path, "create_rmd", "{target}_{ref}_{treat}_differential_signal.log")
    threads: 1
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_differential_rmd.R"

