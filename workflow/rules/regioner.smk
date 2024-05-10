rule localz_regions:
    input:
        checks = ALL_CHECKS,
        features = os.path.join(annotation_path, "features.rds"),
        peaks = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
        ),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
    output:
        rds = os.path.join(peak_path, "{target}", "{target}_regions_localz.rds")
    params:
        regioner_params = extra_params['regioner']           
    threads: 8
    retries: 1
    resources:
        mem_mb = 32768,
        run_time = "30m",
    log: os.path.join(log_path, "regioner", "{target}_regions_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_regions.R"

rule shared_localz_regions:
    input:
        checks = ALL_CHECKS,
        features = os.path.join(annotation_path, "features.rds"),
        peaks = os.path.join(peak_path, "shared", "shared_peaks.bed.gz"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
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
        checks = ALL_CHECKS,
        peaks = expand(
            os.path.join(
                peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
            ),
            target = targets
        ),
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_targets.R"
        ),
        sq = os.path.join(annotation_path, "seqinfo.rds")
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

rule localz_regions_dsa:
    input:
        checks = ALL_CHECKS,
        features = os.path.join(annotation_path, "features.rds"),
        peaks = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-changed.bed.gz"
        ),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_regions.R"
        ),
    output:
        rds = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-regions_localz.rds"
        )
    params:
        regioner_params = extra_params['regioner']        
    threads: 8
    retries: 1
    resources:
        mem_mb = 32768,
        run_time = "30m",
    log: os.path.join(log_path, "regioner", "{target}_{ref}_{treat}_regions_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_regions.R"

rule localz_regions_pairwise:
    input:
        checks = ALL_CHECKS,
        features = os.path.join(annotation_path, "features.rds"),
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise-results.rds"
        ),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_pairwise.R"
        ),
    output:
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise_localz.rds"
        )
    params:
        regioner_params = extra_params['regioner']
    threads: 16
    retries: 1
    resources:
        mem_mb = 64000,
        run_time = "2h",
    log: os.path.join(log_path, "regioner", "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_pairwise.R"