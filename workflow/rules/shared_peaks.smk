rule make_shared_consensus_peaks:
    input:
        arg_checks = rules.check_args.output,
        features = rules.prep_features.output.rds, 
        gtf = rules.create_genome_annotations.output.gtf,
        hic = rules.prep_hic.output.hic, 
        regions = rules.create_genome_annotations.output.regions, 
        peaks = expand(
            os.path.join(
                peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
            ),
            target = targets
        ),
        script = os.path.join("workflow", "scripts", "make_consensus_peaks.R"),
        sq = rules.create_genome_annotations.output.seqinfo, 
        yaml = os.path.join("config", "params.yml"),
    output:
        bed = os.path.join(peak_path, "shared", "shared_peaks.bed.gz"),
        rds =  os.path.join(peak_path, "shared", "shared_peaks.rds"),    
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1
    log: os.path.join(log_path, "make_consensus_peaks", "shared.log")
    resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_shared_consensus_peaks.R"

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

rule localz_regions_shared:
    input:
        arg_checks = rules.check_args.output,
        features = rules.prep_features.output.rds, 
        packages = rules.check_r_packages.output,
        peaks = os.path.join(peak_path, "shared", "shared_peaks.bed.gz"),
        regions = rules.create_genome_annotations.output.regions, 
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_regions.R"
        ),
    output:
        rds = os.path.join(peak_path, "shared", "shared_regions_localz.rds")
    params:
        regioner_params = extra_params['regioner']           
    threads: 8
    retries: 1
    resources:
        mem_mb = 32768,
        run_time = "60m",
    log: os.path.join(log_path, "regioner", "shared_regions_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_regions.R"

rule localz_targets:
    input:
        arg_checks = rules.check_args.output,
        packages = rules.check_r_packages.output,
        peaks = expand(
            os.path.join(
                peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
            ),
            target = targets
        ),
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_targets.R"
        ),
        sq = rules.create_genome_annotations.output.seqinfo, 
    output:
        rds = os.path.join(peak_path, "shared", "shared_targets_localz.rds")
    params:
        regioner_params = extra_params['regioner']           
    threads: 16
    retries: 1
    resources:
        mem_mb = 65536,
        run_time = "2h",
    log: os.path.join(log_path, "regioner", "shared_targets_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_targets.R"