
rule make_consensus_peaks:
    input:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        greylist = lambda wildcards: expand(
            os.path.join(annotation_path, "{f}_greylist.bed.gz"),
            f = set(df[df.target == wildcards.target]['input'])
        ),
        peaks = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}",
                "{{target}}_{treat}_filtered_peaks.narrowPeak"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
    output:
        peaks = os.path.join(
            macs2_path, "{target}", "{target}_consensus_peaks.bed.gz"
        )
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
        peaks = expand(
            os.path.join(
                macs2_path, "{target}", "{target}_consensus_peaks.bed.gz"
            ),
            target = [targets]
        ),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
    output:
        peaks = os.path.join(
            macs2_path, "shared", "shared_consensus_peaks.bed.gz"
        )
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "make_consensus_peaks", "shared.log")
        resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_shared_consensus_peaks.R"