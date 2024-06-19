rule count_windows:
    input:
        arg_checks = rules.check_args.output,
        bam = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.bam"),
            sample = df['sample'][(df['target'] == wildcards.target)]
        ),
        bai = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.bam.bai"),
            sample = df['sample'][(df['target'] == wildcards.target)]
        ),
        input_bam = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.bam"),
            sample = set(df['input'][(df['target'] == wildcards.target)])
        ),
        input_bai = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.bam.bai"),
            sample = set(df['input'][(df['target'] == wildcards.target)])
        ),
        blacklist = rules.prep_blacklist.output.blacklist, 
        greylist = rules.combine_greylists.output.rds,
        macs2_logs = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}",
                "{{target}}_{treat_levels}_merged_callpeak.log"
            ),
            treat_levels = set(df['treat'][df['target'] == wildcards.target])
        ),
        packages = rules.check_r_packages.output,   
        peak_qc = os.path.join(
            macs2_path, "{target}", "{target}_qc_samples.tsv"
        ),
        peaks = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.rds"
        ),
        script = os.path.join("workflow", "scripts", "make_counts.R"),
        seqinfo = rules.create_genome_annotations.output.seqinfo, 
    output:
        rds = os.path.join(diff_path, "{target}", "{target}_counts.rds")
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "count_windows", "{target}_make_counts.log")
    retries: 1
    threads: lambda wildcards, attempt: attempt * 8
    params:
        contrasts = lambda wildcards: diff_sig_param[wildcards.target]['contrasts'],
        filter_q = lambda wildcards: diff_sig_param[wildcards.target]['filter_q'],
        win_type = lambda wildcards: diff_sig_param[wildcards.target]['window_type'],
        win_size = lambda wildcards: diff_sig_param[wildcards.target]['window_size'],
        win_step = lambda wildcards: diff_sig_param[wildcards.target]['window_step'],
    resources:
        runtime = lambda wildcards, attempt: attempt * 60,
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    script:
        "../scripts/make_counts.R"

rule differential_signal_analysis:
    input:
        arg_checks = rules.check_args.output,
        counts = os.path.join(diff_path, "{target}", "{target}_counts.rds"),
        gtf = rules.create_genome_annotations.output.gtf,
        hic = rules.prep_hic.output.hic, 
        features = rules.prep_features.output.rds, 
        packages = rules.check_r_packages.output,   
        peaks = expand(
            os.path.join(
                peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
                ),
                target = targets
        ),
        regions = rules.create_genome_annotations.output.regions, 
        script = os.path.join("workflow", "scripts", "differential_signal.R"),
        sq = rules.create_genome_annotations.output.seqinfo, 
        yaml = os.path.join("config", "params.yml"),
    output:
        changed = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-changed.bed.gz"
        ),
        decreased = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-decreased.bed.gz"
        ),
        increased = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-increased.bed.gz"
        ),
        ihw = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-ihw.rds"
        ),
        rds = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}-differential-signal.rds"
        ),
    params:
        diff_sig_params = lambda wildcards: diff_sig_param[wildcards.target],
        peak_calling_params = lambda wildcards: peak_calling_param[wildcards.target],
    threads: 6
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "differential_signal", "{target}_{ref}_{treat}.log")
    resources:
        runtime = "30m",
        mem_mb = 64000,
    script:
        "../scripts/differential_signal.R"

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

