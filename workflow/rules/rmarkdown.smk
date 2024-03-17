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

rule compile_index_html:
    input:
        html = HTML_OUT,
        rmd = os.path.join(rmd_path, "index.Rmd"),
        setup = rules.create_setup_chunk.output,
        site_yaml = rules.create_site_yaml.output,
        rulegraph = 'workflow/rules/rulegraph.dot'
    output:
        html = "docs/index.html"
    conda: "../envs/rmarkdown.yml"
    threads: 1
    resources:
        mem_mb = 1024,
        runtime = "5m",
    log: os.path.join(log_path, "rmarkdown/compile_index_html.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule create_macs2_summary_rmd:
    input:
        here = rules.check_here_file.output,
        packages = rules.check_r_packages.output,
        module = "workflow/modules/macs2_summary.Rmd",
    output:
        rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd")
    params:
        min_prop = lambda wildcards: macs2_qc_param[wildcards.target]['min_prop_reps'],
        outlier_thresh = lambda wildcards: macs2_qc_param[wildcards.target]['outlier_threshold'],
        macs2_fdr = lambda wildcards: macs2_param[wildcards.target]['fdr'],
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "macs2_summary", "create_{target}_macs2_summary.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_macs2_summary.R"


rule compile_macs2_summary_html:
    input:
        annotations = ALL_RDS,
        assets = os.path.join("docs", "assets"),
        bw = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}",
                "{{target}}_{treat}_merged_treat_pileup.bw"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        external = rules.check_external_files.output,
        cors = os.path.join(
            macs2_path, "{target}", "{target}_cross_correlations.tsv"
        ),
        peak_files = expand(
             os.path.join(
                peak_path, "{{target}}", "{{target}}_{f}"
            ),
            f = [
                'motif_position.tsv.gz', 'motif_enrichment.tsv.gz', 
                'regions_localz.rds', 'consensus_peaks.bed.gz'
                ]
        ),
        rmd = os.path.join(rmd_path, "{target}_macs2_summary.Rmd"),
        setup = rules.create_setup_chunk.output,
        yaml = rules.create_site_yaml.output
    output:
        html = os.path.join("docs", "{target}_macs2_summary.html"),
        fig_path = directory(
            os.path.join("docs", "{target}_macs2_summary_files", "figure-html")
        ),
        great = os.path.join(
            "output", "results", "{target}", "{target}_great_results.tsv.gz"
        ),
        localz = os.path.join(
            "output", "results", "{target}", "{target}_localz.tsv"
        ),
        renv = temp(
            os.path.join("output", "envs", "{target}_macs2_summary.RData")
        ),
    conda: "../envs/rmarkdown.yml"
    threads: 6
    resources:
        mem_mb = 16384,
        runtime = "30m",
    log: os.path.join(log_path, "macs2_summary", "compile_{target}_macs2_summary.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule compile_peak_comparison_rmd:
    input:
        annotations = ALL_RDS,
        motif_enrich = os.path.join(
            peak_path, "shared", "shared_motif_enrichment.tsv.gz"
        ),
        motif_pos = os.path.join(
            peak_path, "shared", "shared_motif_position.tsv.gz"
        ),
        rmd = os.path.join("workflow", "modules", "peak_comparison.Rmd"),
        shared_files = expand(
            os.path.join(peak_path, "shared", "shared_{f}"),
            f = [
                'consensus_peaks.bed.gz', 'motif_position.tsv.gz', 
                'motif_enrichment.tsv.gz', 'targets_localz.rds',
                'regions_localz.rds',
                ]
        ),
        target_files = expand(
            os.path.join(peak_path, "{tg}", "{tg}_{f}"),
            tg = targets,
            f = ['regions_localz.rds', 'consensus_peaks.bed.gz']
        ),
    output:
        rmd = os.path.join(rmd_path, "peak_comparison.Rmd"),
        html = os.path.join("docs", "peak_comparison.html"),
        tsv = expand(
            os.path.join("output", "results", "shared", "{f}"),
            f = ['shared_great_results_genomic_bg.tsv.gz',
            'shared_great_results_targets_bg.tsv.gz',
            'shared_regions_localz.tsv', 'pairwise_localz.tsv']
        )
    conda: "../envs/rmarkdown.yml"
    threads: 6
    resources:
        mem_mb = 16384,
        runtime = "30m",
    log: os.path.join(log_path, "peak_comparison", "compile_peak_comparison.log")
    shell:
        """
        cp {input.rmd} {output.rmd}
        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}
        """
