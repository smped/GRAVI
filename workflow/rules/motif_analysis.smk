rule motif_analysis_peaks:
    input:
        arg_checks = rules.check_args.output,
        exclude_ranges = rules.make_exclude_ranges.output.rds,
        gene_regions = rules.create_genome_annotations.output.regions, 
        motifs = rules.prep_motifs.output.motifs,
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

rule motif_analysis_shared:
    input:
        arg_checks = rules.check_args.output,
        exclude_ranges = rules.make_exclude_ranges.output.rds,
        gene_regions = rules.create_genome_annotations.output.regions, 
        motifs = rules.prep_motifs.output.motifs,
        packages = rules.check_r_packages.output,
        peaks = os.path.join(peak_path, "shared", "shared_peaks.rds"),
        script = os.path.join("workflow", "scripts", "motif_analysis.R"),
    output:
        enrich = os.path.join(
            peak_path, "shared", "shared_motif_enrichment.tsv.gz"
        ),
        pos = os.path.join(
            peak_path, "shared", "shared_motif_position.tsv.gz"
        ),
    params:
        motif_params = motif_param['shared'],
    threads: lambda wildcards, attempt: attempt * 8
    retries: 1
    resources:
        disk_mb = 10000,
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        runtime = lambda wildcards, attempt: attempt * 120,
    log: os.path.join(log_path, "motif_analysis", "shared_motif_analysis.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/motif_analysis.R"

rule motif_analysis_dsa:
    input:
        arg_checks = rules.check_args.output,
        motifs = rules.prep_motifs.output.motifs,
        packages = rules.check_r_packages.output,
        results = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}-differential-signal.rds"
        ),
        script = os.path.join("workflow", "scripts", "motif_analysis_dsa.R"),
    output:
        enrich = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}_motif_enrichment.tsv.gz"
        ),
        pos = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}_motif_position.tsv.gz"
        ),
    params:
        motif_params = lambda wildcards: motif_param[wildcards.target],
    threads: lambda wildcards, attempt: attempt * 4
    retries: 2
    resources:
        disk_mb = 10000,
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        runtime = lambda wildcards, attempt: attempt * 30,
    log: os.path.join(log_path, "motif_analysis_dsa", "{target}_{ref}_{treat}.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/motif_analysis_dsa.R"

rule motif_analysis_pairwise:
    input:
        arg_checks = rules.check_args.output,
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise-results.rds"
        ),
        motifs = rules.prep_motifs.output.motifs,
        packages = rules.check_r_packages.output,
        script = os.path.join(
            "workflow", "scripts", "pairwise_motif_analysis.R"
        ),
        seqinfo = rules.create_genome_annotations.output.seqinfo, 
    output:
        enrich_tsv = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-motif_enrichment.tsv.gz"
        ),
        position_tsv = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-motif_position.tsv.gz"
        ),        
    params:
        motif_params = motif_param['pairwise']
    retries: 1
    threads: lambda wildcards, attempt: attempt * 8
    conda: "../envs/rmarkdown.yml"
    resources:
        disk_mb = 10000,
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        runtime = lambda wildcards, attempt: attempt * 30,
    log: os.path.join(log_path, "pairwise_motif_analysis", "{tgt1}_{comp1}-{tgt2}_{comp2}.log")
    script:
        "../scripts/pairwise_motif_analysis.R"    

