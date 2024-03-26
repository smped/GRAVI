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
    log: os.path.join(log_path, "compile_rmd", "compile_index_html.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """

rule compile_annotations_html:
    input:
        checks = ALL_CHECKS,
        greylist = expand(
            os.path.join(grey_path, "{f}_greylist.bed.gz"),
            f = set(df['input'])
        ),
        rmd = os.path.join(rmd_path, "annotation_description.Rmd"),
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
    log: os.path.join(log_path, "compile_rmd", "compile_annotations_html.log")
    resources:
        mem_mb = 4096,
        disk_mb = 4000,
        run_time = "10m",
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """


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
    log: os.path.join(log_path, "compile_rmd", "compile_{target}_macs2_summary.log")
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
    output:
        rmd = os.path.join(rmd_path, "peak_comparison.Rmd"),
        html = os.path.join("docs", "peak_comparison.html"),
        tsv = expand(
            os.path.join("output", "results", "shared", "{f}"),
            f = ['shared_enrichment_results_genomic_bg.tsv.gz',
            'shared_enrichment_results_targets_bg.tsv.gz',
            'shared_regions_localz.tsv', 'pairwise_localz.tsv']
        )
    conda: "../envs/rmarkdown.yml"
    threads: 6
    resources:
        mem_mb = 16384,
        runtime = "30m",
    log: os.path.join(log_path, "compile_rmd", "compile_peak_comparison.log")
    shell:
        """
        cp {input.rmd} {output.rmd}
        R -e "rmarkdown::render_site('{output.rmd}')" &>> {log}
        """

rule compile_differential_signal_html:
    input:
        annotations = ALL_RDS,
        counts = os.path.join(diff_path, "{target}", "{target}_counts.rds"),
        rmd = os.path.join(
            rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
        ),
    output:
        html = "docs/{target}_{ref}_{treat}_differential_signal.html",
        enrichment = expand(
            os.path.join(
                diff_path, "{{target}}",
                "{{target}}_{{ref}}_{{treat}}-{f}-enrichment.csv"
            ),
            f = ['changed', 'increased', 'decreased']
        ),
        fig_path = directory(
            os.path.join(
                "docs", "{target}_{ref}_{treat}_differential_signal_files",
                "figure-html"
            )
        ),
        results = os.path.join(
            diff_path, "{{target}}",
            "{{target}}_{{ref}}_{{treat}}-differential_signal.csv.gz"
        ),
        renv = temp(
            os.path.join(
                "output", "envs",
                "{target}_{ref}_{treat}-differential_signal.RData"
            )
        ),
    conda: "../envs/rmarkdown.yml"
    threads: 8
    resources:
        mem_mb = 65536,
        runtime = "1h"
    log: os.path.join(log_path, "compile_rmd", "{target}_{ref}_{treat}_differential_signal.log")
    shell:
        """
        R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
        """
