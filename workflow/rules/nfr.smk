rule merge_filtered_peaks:
    input:
        peaks = os.path.join(
            peak_path, "{target}", "{target}_{treat}_filtered_peaks.narrowPeak"
        ),
        script = os.path.join("workflow", "scripts", "merge_filtered_peaks.R"),
        yaml = os.path.join("config", "params.yml")
    output:
        bed = os.path.join(
            peak_path, "{target}", "{target}_{treat}_merged_peaks.bed"
        )
    params:
        within = 300
    threads: 1
    conda: "../envs/rmarkdown.yml"
    resources:
        mem_mb = 8192,
    log: os.path.join(log_path, "merge_filtered_peaks", "{target}_{treat}.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/merge_filtered_peaks.R"

rule call_nfr:
    input:
        bdg = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_treat_pileup.bdg"
        ),
        bed =  os.path.join(
            peak_path, "{target}", "{target}_{treat}_merged_peaks.bed"
        ),
        script = os.path.join("workflow", "scripts", "HisTrader.pl")
    output:
        nfr = temp(
            os.path.join(
                peak_path, "{target}", "{target}_{treat}.nfr.bed"
            )
        ),
    params:
        pre = "{target}_{treat}",
        p_max = 0.1,
        max_nfr = 1000,
        min_size = 500,
    threads: 1
    resources:
        mem_mb = 8192,
        runtime = "1h"
    shadow: 'minimal'
    log: os.path.join(log_path, "call_nfr", "{target}_{treat}.log")
    shell:
        """
        perl {input.script} \
          --bedGraph {input.bdg} \
          --peaks {input.bed} \
          --pMax {params.p_max} \
          --minSize {params.min_size} \
          --filter {params.max_nfr} \
          --out {params.pre}
        """

rule strip_nfr_bed:
    input: "{f}.bed"
    output: "{f}.bed.gz"
    threads: 1
    resources:
        runtime = "5m",
        mem_mb = 2048
    shell:
        """
        cut -f1-3 {input} | gzip -c > {output}
        """

rule make_consensus_nfr:
    input:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        features = os.path.join(annotation_path, "features.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        greylist = lambda wildcards: expand(
            os.path.join(grey_path, "{f}_greylist.bed.gz"),
            f = set(df[df.target == wildcards.target]['input'])
        ),
        hic = os.path.join(annotation_path, "hic.rds"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join("workflow", "scripts", "make_consensus_peaks.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
        yaml = os.path.join("config", "params.yml"),
        peaks = lambda wildcards: expand(
            os.path.join(
                peak_path, "{{target}}", "{{target}}_{treat}.nfr.bed.gz"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        )
    output:
        bed = os.path.join(
            peak_path, "{target}", "{target}_consensus_nfr.bed.gz"
        ),
        rds = os.path.join(
            peak_path, "{target}", "{target}_consensus_nfr.rds"
        )
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "make_consensus_nfr", "{target}.log")
    resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_consensus_peaks.R"
    