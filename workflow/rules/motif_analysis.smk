rule motif_analysis_peaks:
    input:
        arg_checks = rules.check_args.output,
        exclude_ranges = rules.make_exclude_ranges.output.rds,
        gene_regions = rules.create_genome_annotations.output.regions, 
        motifs = rules.prep_motifs.output.motifs,
        n_masked_ranges = rules.make_n_masked_ranges.output.rds,
        packages = rules.check_r_packages.output,
        peaks = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.rds"
        ),
        script = os.path.join("workflow", "scripts", "motif_analysis.R"),
    output:
        enrich = os.path.join(
            peak_path, "{target}", "{target}_motif_enrichment.tsv.gz"
        ),
        pos = os.path.join(
            peak_path, "{target}", "{target}_motif_position.tsv.gz"
        ),
    params:
        motif_params = lambda wildcards: motif_param[wildcards.target]
    threads: lambda wildcards, attempt: attempt * 8
    retries: 2
    resources:
        disk_mb = 10000,
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        runtime = lambda wildcards, attempt: attempt * 120,
    log: os.path.join(log_path, "motif_analysis", "{target}_motif_analysis.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/motif_analysis.R"


