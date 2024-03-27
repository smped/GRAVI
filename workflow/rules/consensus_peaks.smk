rule filter_merged_peaks:
    input:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        greylist = lambda wildcards: expand(
            os.path.join(grey_path, "{f}_greylist.bed.gz"),
            f = set(df[df.target == wildcards.target]['input'])
        ),
        merged = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_peaks.narrowPeak"
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv"),
        rep = lambda wildcards: expand(
            os.path.join(macs2_path, "{f}", "{f}_peaks.narrowPeak"),
            f = set(df[(df.treat == wildcards.treat) & (df.target == wildcards.target)]['sample'])
        ),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
    output:
        peaks = os.path.join(
            peak_path, "{target}", "{target}_{treat}_filtered_peaks.narrowPeak"
        )
    params:
        min_prop = lambda wildcards: peak_qc_param[wildcards.target]['min_prop_reps']
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "filter_merged_peaks", "{target}_{treat}.log")
    resources:
        mem_mb = 4096,
        runtime = "15m"
    script:
        "../scripts/filter_merged_peaks.R"

rule make_consensus_peaks:
    input:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        features = os.path.join(annotation_path, "features.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        greylist = lambda wildcards: expand(
            os.path.join(grey_path, "{f}_greylist.bed.gz"),
            f = set(df[df.target == wildcards.target]['input'])
        ),
        hic = os.path.join(annotation_path, "hic.rds"),
        peaks = lambda wildcards: expand(
            os.path.join(
                peak_path, "{{target}}",
                "{{target}}_{treat}_filtered_peaks.narrowPeak"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join("workflow", "scripts", "make_consensus_peaks.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
        yaml = os.path.join("config", "params.yml"),
    output:
        peaks = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
        ),
        rds =  os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.rds"
        ),
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "make_consensus_peaks", "{target}.log")
    resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_consensus_peaks.R"

rule make_shared_consensus_peaks:
    input:
        features = os.path.join(annotation_path, "features.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        hic = os.path.join(annotation_path, "hic.rds"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        peaks = expand(
            os.path.join(
                peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
            ),
            target = targets
        ),
        script = os.path.join("workflow", "scripts", "make_consensus_peaks.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
        yaml = os.path.join("config", "params.yml"),
    output:
        peaks = os.path.join(peak_path, "shared", "shared_peaks.bed.gz"),
        rds =  os.path.join(peak_path, "shared", "shared_peaks.rds"),        
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "make_consensus_peaks", "shared.log")
    resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_shared_consensus_peaks.R"