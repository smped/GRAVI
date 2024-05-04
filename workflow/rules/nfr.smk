rule merge_nearby_peaks:
    input:
        peaks = os.path.join(
            peak_path, "{target}", "{target}_{treat}_filtered_peaks.narrowPeak"
        ),
        script = os.path.join("workflow", "scripts", "merge_filtered_peaks.R"),
    output:
        bed = os.path.join(
            nfr_path, "{target}", "{target}_{treat}_merged_peaks.bed"
        )
    params:
        within = nfr_params['merge_peaks_within']
    threads: 1
    conda: "../envs/rmarkdown.yml"
    resources:
        mem_mb = 8192,
        runtime = "5m"
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
            nfr_path, "{target}", "{target}_{treat}_merged_peaks.bed"
        ),
        script = os.path.join("workflow", "scripts", "HisTrader.pl")
    output:
        nfr = temp(
            os.path.join(nfr_path, "{target}", "{target}_{treat}.nfr.bed")
        ),
    params:
        pre = os.path.join(nfr_path, "{target}", "{target}_{treat}"),
        p_max = nfr_params['p_max'],
        max_nfr = max(nfr_params['nfr_width']),
        min_size = nfr_params['min_peak_width'],
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
          --out {params.pre} > {log}
        """

rule strip_bed:
    input: "{f}.bed"
    output: "{f}.bed.gz"
    threads: 1
    localrule: True
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
        greylist = os.path.join(grey_path, "greylists.rds"),
        hic = os.path.join(annotation_path, "hic.rds"),
        peaks = lambda wildcards: expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_{treat}.nfr.bed.gz"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join("workflow", "scripts", "make_consensus_peaks.R"),
        sq = os.path.join(annotation_path, "seqinfo.rds"),
        yaml = os.path.join("config", "params.yml"),
    output:
        bed = os.path.join(
            nfr_path, "{target}", "{target}_consensus_nfr.bed.gz"
        ),
        rds = os.path.join(
            nfr_path, "{target}", "{target}_consensus_nfr.rds"
        )
    params:
        ## Passed to makeConsensus. This will give stringent, shared NFRs
        method = 'coverage',
        min_width = min(nfr_params['nfr_width']),
        p = 1,
        min_gapwidth = nfr_params['merge_nfrs_within']
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "make_consensus_peaks", "{target}_nfr.log")
    resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_consensus_peaks.R"

rule nfr_motif_analysis:
    input:
        exclude_ranges = os.path.join(annotation_path, "exclude_ranges.rds"),
        gene_regions = os.path.join(annotation_path, "gene_regions.rds"),
        motifs = os.path.join(annotation_path, "motif_list.rds"),
        packages = os.path.join(check_path, "r-packages.chk"),
        peaks = os.path.join(
            nfr_path, "{target}", "{target}_consensus_nfr.rds"
        ),
        script = os.path.join("workflow", "scripts", "motif_analysis.R"),
    output:
        enrich = os.path.join(
            nfr_path, "{target}", "{target}_motif_enrichment.tsv.gz"
        ),
        pos = os.path.join(
            nfr_path, "{target}", "{target}_motif_position.tsv.gz"
        ),
    params:
        motif_params = motif_param['nfr']
    threads: lambda wildcards, attempt: attempt * 8
    retries: 2
    resources:
        disk_mb = 10000,
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        runtime = lambda wildcards, attempt: attempt * 120,
    log: os.path.join(log_path, "motif_analysis", "{target}_nfr_motif_analysis.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/motif_analysis.R"

rule nfr_localz_regions:
    input:
        checks = ALL_CHECKS,
        features = os.path.join(annotation_path, "features.rds"),
        peaks = os.path.join(
            nfr_path, "{target}", "{target}_consensus_nfr.bed.gz"
        ),
        params = os.path.join("config", "params.yml"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_regions.R"
        ),
    output:
        rds = os.path.join(
            nfr_path, "{target}", "{target}_nfr_regions_localz.rds"
        )
    threads: 8
    retries: 1
    resources:
        mem_mb = 32768,
        run_time = "60m",
    log: os.path.join(log_path, "regioner", "{target}_nfr_regions_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_regions.R"

def get_nfr_peaks_for_local_z(wildcards):
    tgts = set(targets).difference(set([wildcards.target]))
    peaks = []
    peaks.extend(
        expand(
            os.path.join(nfr_path, "{t}", "{t}_consensus_nfr.bed.gz"),
            t = [wildcards.target]
        )
    )
    peaks.extend(
        expand(
            os.path.join(peak_path, "{t}", "{t}_consensus_peaks.bed.gz"),
            t = set(targets).difference(set([wildcards.target]))
        )
    )
    return(peaks)


rule nfr_localz_targets:
    input:
        checks = ALL_CHECKS,
        params = os.path.join("config", "params.yml"),
        peaks = get_nfr_peaks_for_local_z,
        script = os.path.join(
            "workflow", "scripts", "regioner_localz_targets.R"
        ),
        sq = os.path.join(annotation_path, "seqinfo.rds")
    output:
        rds = os.path.join(
            nfr_path, "{target}", "{target}_nfr_targets_localz.rds"
        )
    threads: 16
    retries: 1
    resources:
        mem_mb = 65536,
        run_time = "2h",
    log: os.path.join(log_path, "regioner", "{target}_nfr_targets_localz.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/regioner_localz_targets.R"
