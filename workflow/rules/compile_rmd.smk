rule compile_index_html:
    input:
        html = HTML_OUT,
        packages = rules.check_r_packages.output,
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
    log: os.path.join(log_path, "compile_rmd", "compile_index_html.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule compile_annotations_html:
    input:
        blacklist = blacklist,
        rmd = os.path.join(rmd_path, "annotation_description.Rmd"),
        rna_module = os.path.join(
            "workflow", "modules", "_rna_description.Rmd"
        ),        
        setup = rules.create_setup_chunk.output,
        site_yaml = rules.create_site_yaml.output
    output:
        rds = os.path.join(annotation_path, "colours.rds"),
        html = "docs/annotation_description.html",
        fig_path = directory(
            os.path.join("docs", "annotation_description_files", "figure-html")
        )
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1
    log: os.path.join(log_path, "compile_rmd", "compile_annotations_html.log")
    resources:
        mem_mb = 4096,
        disk_mb = 4000,
        run_time = "10m",
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """


rule compile_signal_summary_html:
    input:
        annotations = ANNOTATION_RDS,
        rmd = os.path.join(rmd_path, "{target}_signal_summary.Rmd"),
        setup = rules.create_setup_chunk.output,
        yaml = rules.create_site_yaml.output
    output:
        html = os.path.join("docs", "{target}_signal_summary.html"),
        fig_path = directory(
            os.path.join("docs", "{target}_signal_summary_files", "figure-html")
        ),
        great = os.path.join(
            "output", "results", "{target}", "{target}_enrichment.tsv"
        ),
        localz = os.path.join(
            "output", "results", "{target}", "{target}_localz.tsv"
        ),
        renv = temp(
            os.path.join("output", "envs", "{target}_signal_summary.RData")
        ),
    conda: "../envs/rmarkdown.yml"
    threads: 6
    retries: 0
    resources:
        mem_mb = 16384,
        runtime = "30m",
    log: os.path.join(log_path, "compile_rmd", "compile_{target}_signal_summary.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule compile_signal_comparison_html:
    input:
        annotations = ANNOTATION_RDS,
        arg_checks = rules.check_args.output,
        bigwig = MERGED_BW,
        motif_enrich = os.path.join(
            peak_path, "shared", "shared_motif_enrichment.tsv.gz"
        ),
        motif_pos = os.path.join(
            peak_path, "shared", "shared_motif_position.tsv.gz"
        ),
        nfr = NFR_RDS,
        packages = rules.check_r_packages.output,
        rmd = os.path.join("workflow", "modules", "signal_comparison.Rmd"),
        shared_files = expand(
            os.path.join(peak_path, "shared", "shared_{f}"),
            f = [
                'peaks.bed.gz', 'motif_position.tsv.gz',
                'motif_enrichment.tsv.gz', 'targets_localz.rds',
                'regions_localz.rds',
                ]
        ),
        target_files = expand(
            os.path.join(peak_path, "{tg}", "{tg}_{f}"),
            tg = targets,
            f = ['regions_localz.rds', 'consensus_peaks.bed.gz']
        ),
        setup = rules.create_setup_chunk.output,
        yaml = rules.create_site_yaml.output
    output:
        rmd = os.path.join(rmd_path, "signal_comparison.Rmd"),
        fig_path = directory(
            os.path.join("docs", "signal_comparison_files", "figure-html")
        ),        
        html = os.path.join("docs", "signal_comparison.html"),
        tsv = expand(
            os.path.join("output", "results", "shared", "{f}"),
            f = ['shared_enrichment_genomic_bg.tsv',
            'shared_enrichment_targets_bg.tsv',
            'shared_regions_localz.tsv', 'pairwise_localz.tsv']
        )
    conda: "../envs/rmarkdown.yml"
    threads: 6
    retries: 1
    resources:
        mem_mb = 16384,
        runtime = "30m",
    log: os.path.join(log_path, "compile_rmd", "signal_comparison.log")
    shell:
        """
        cp {input.rmd} {output.rmd}
        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}
        """

rule compile_nfr_html:
    input:
        annotations = ANNOTATION_RDS,
        rmd = os.path.join(rmd_path, "{target}_nfr.Rmd"),
        setup = rules.create_setup_chunk.output,
        yaml = rules.create_site_yaml.output
    output:
        html = os.path.join("docs", "{target}_nfr.html"),
        fig_path = directory(
            os.path.join("docs", "{target}_nfr_files", "figure-html")
        ),
    conda: "../envs/rmarkdown.yml"
    retries: 1
    threads: lambda wildcards, attempt: 8 * attempt
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        runtime = "30m",
    log: os.path.join(log_path, "compile_rmd", "{target}_nfr.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule compile_differential_signal_html:
    input:
        modules = expand(
            os.path.join("workflow", "modules", "_{f}_differential_signal.Rmd"),
            f = ['rna', 'ihw', 'nfr']
        ),
        rmd = os.path.join(
            rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
        ),
        setup = rules.create_setup_chunk.output,
        yaml = rules.create_site_yaml.output
    output:
        asset_path = directory(
            os.path.join("docs", "assets", "{target}_{ref}_{treat}")
        ),
        html = "docs/{target}_{ref}_{treat}_differential_signal.html",
        # enrichment = expand(
        #     os.path.join(
        #         diff_path, "{{target}}",
        #         "{{target}}_{{ref}}_{{treat}}-{f}-enrichment.csv"
        #     ),
        #     f = ['changed', 'increased', 'decreased']
        # ),
        fig_path = directory(
            os.path.join(
                "docs", "{target}_{ref}_{treat}_differential_signal_files",
                "figure-html"
            )
        ),
        # results = os.path.join(
        #     diff_path, "{{target}}",
        #     "{{target}}_{{ref}}_{{treat}}-differential_signal.csv.gz"
        # ),
        renv = temp(
            os.path.join(
                "output", "envs",
                "{target}_{ref}_{treat}-differential_signal.RData"
            )
        ),
    conda: "../envs/rmarkdown.yml"
    retries: 1
    threads: lambda wildcards, attempt: 6 * attempt
    resources:
        mem_mb = lambda wildcards, attempt: 48000 * attempt,
        runtime = lambda wildcards, attempt: 30 * attempt,
        disk_mb = lambda wildcards, attempt: 5000 * attempt,
    log: os.path.join(log_path, "compile_rmd", "{target}_{ref}_{treat}_differential_signal.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule compile_pairwise_comparison_html:
    input:
        annotations = ANNOTATION_RDS,
        results = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise-results.rds"
        ),
        rmd = os.path.join(
            rmd_path, "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise_comparison.Rmd"
        ),
    output:
        html = os.path.join(
            "docs", "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise_comparison.html"
        ),
        fig_path = directory(
            os.path.join(
                "docs", 
                "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise_comparison_files",
                "figure-html"
            )
        )
    conda: "../envs/rmarkdown.yml"
    threads: 4
    resources:
        mem_mb = 32000,
        runtime = "30m"
    log: os.path.join(log_path, "compile_rmd", "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise_comparison.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """