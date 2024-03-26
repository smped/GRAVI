rule count_windows:
    input:
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
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        chk = ALL_CHECKS,
        greylist = lambda wildcards: expand(
            os.path.join(grey_path, "{ip_sample}_greylist.bed.gz"),
            ip_sample = set(df['input'][df['target'] == wildcards.target])
        ),
        macs2_logs = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}",
                "{{target}}_{treat_levels}_merged_callpeak.log"
            ),
            treat_levels = set(df['treat'][df['target'] == wildcards.target])
        ),
        macs2_qc = os.path.join(
            macs2_path, "{target}", "{target}_qc_samples.tsv"
        ),
        peaks = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.rds"
        ),
        script = os.path.join("workflow", "scripts", "make_counts.R"),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
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
        counts = os.path.join(diff_path, "{target}", "{target}_counts.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        hic = os.path.join(annotation_path, "hic.rds"),
        features = os.path.join(annotation_path, "features.rds"),
        peaks = expand(
            os.path.join(
                peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
                ),
                target = targets
        ),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join("workflow", "scripts", "differential_signal.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
        yaml = os.path.join("config", "params.yml"),
    output:
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
        diff_sig_params = lambda wildcards: diff_sig_param[wildcards.target]
    threads: 6
    conda: "../envs/rmarkdown.yml"
    log: os.path.join(log_path, "differential_signal", "{target}_{ref}_{treat}.log")
    resources:
        runtime = "1h",
        mem_mb = 64000,
    script:
        "../scripts/differential_signal.R"
