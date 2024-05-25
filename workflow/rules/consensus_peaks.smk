rule filter_merged_peaks:
    input:
        blacklist = rules.prep_blacklist.output.blacklist, 
        greylist = rules.combine_greylists.output.rds,
        merged = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_peaks.narrowPeak"
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv"),
        rep = lambda wildcards: expand(
            os.path.join(macs2_path, "{f}", "{f}_peaks.narrowPeak"),
            f = set(df[(df.treat == wildcards.treat) & (df.target == wildcards.target)]['sample'])
        ),
        sq = rules.create_genome_annotations.output.seqinfo, 
    output:
        peaks = os.path.join(
            peak_path, "{target}", "{target}_{treat}_filtered_peaks.narrowPeak"
        )
    params:
        min_prop = lambda wildcards: peak_calling_param[wildcards.target]['min_prop_reps'],
        merge_fdr = lambda wildcards: peak_calling_param[wildcards.target]['merge_fdr'],
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1
    log: os.path.join(log_path, "filter_merged_peaks", "{target}_{treat}.log")
    resources:
        mem_mb = 4096,
        runtime = "15m"
    script:
        "../scripts/filter_merged_peaks.R"

rule make_consensus_peaks:
    input:
        blacklist = rules.prep_blacklist.output.blacklist, 
        features = rules.prep_features.output.rds, 
        gtf = rules.create_genome_annotations.output.gtf,
        greylist = rules.combine_greylists.output.rds,
        hic = rules.prep_hic.output.hic, 
        peaks = lambda wildcards: expand(
            os.path.join(
                peak_path, "{{target}}",
                "{{target}}_{treat}_filtered_peaks.narrowPeak"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv"),
        regions = rules.create_genome_annotations.output.regions, 
        script = os.path.join("workflow", "scripts", "make_consensus_peaks.R"),
        sq = rules.create_genome_annotations.output.seqinfo, 
        yaml = os.path.join("config", "params.yml"),
    output:
        bed = os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.bed.gz"
        ),
        rds =  os.path.join(
            peak_path, "{target}", "{target}_consensus_peaks.rds"
        ),
    params:
        method = "union",
        merge_within =  lambda wildcards: peak_calling_param[wildcards.target]['merge_within'],
        min_width = lambda wildcards: peak_calling_param[wildcards.target]['min_width'],
        peak_type = lambda wildcards: peak_calling_param[wildcards.target]['peak_type'],
        p = 0,
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1
    log: os.path.join(log_path, "make_consensus_peaks", "{target}.log")
    resources:
        mem_mb = 4096,
        runtime = "10m"
    script:
        "../scripts/make_consensus_peaks.R"

rule make_shared_consensus_peaks:
    input:
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