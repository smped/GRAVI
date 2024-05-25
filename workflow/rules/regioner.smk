rule localz_regions:
    input:
        checks = ALL_CHECKS,
        features = rules.prep_features.output.rds, 
        peaks = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
        ),
        regions = rules.create_genome_annotations.output.regions, 
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

rule localz_regions_shared:
    input:
        checks = ALL_CHECKS,
        features = rules.prep_features.output.rds, 
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

rule localz_regions_dsa:
    input:
        checks = ALL_CHECKS,
        features = rules.prep_features.output.rds, 
        peaks = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-changed.bed.gz"
        ),
        regions = rules.create_genome_annotations.output.regions, 
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
        mem_mb = 32000,
        run_time = "30m",
    log: os.path.join(log_path, "regioner", "{target}_{ref}_{treat}_regions_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_regions.R"

rule localz_regions_pairwise:
    input:
        bed = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-{pw_dir}.bed.gz"
        ),
        checks = ALL_CHECKS,
        features = rules.prep_features.output.rds, 
        regions = rules.create_genome_annotations.output.regions, 
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_pairwise.R"
        ),
    output:
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-{pw_dir}_localz.rds"
        )
    params:
        regioner_params = extra_params['regioner']
    threads: 8
    retries: 1
    resources:
        mem_mb = 32000,
        run_time = "1h",
    log: os.path.join(log_path, "regioner_pairwise", "{tgt1}_{comp1}-{tgt2}_{comp2}_{pw_dir}_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_pairwise.R"

rule merge_localz_pairwise:
    input:
        rds = expand(
            os.path.join(
                pairs_path, "{{tgt1}}_{{comp1}}-{{tgt2}}_{{comp2}}",
                "{{tgt1}}_{{comp1}}-{{tgt2}}_{{comp2}}-{f}_localz.rds"
            ),
            f = pw_dirs
        ),
        script = os.path.join("workflow", "scripts", "merge_localz_pairwise.R")
    output:
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise_localz.rds"
        ),
    threads: 4
    retries: 1
    resources:
        mem_mb = 32000,
        run_time = "1h",
    log: os.path.join(log_path, "regioner_pairwise", "{tgt1}_{comp1}-{tgt2}_{comp2}_merge_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/merge_localz_pairwise.R"        
    
