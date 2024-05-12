rule prepare_pairwise_results:
    input:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        features = os.path.join(annotation_path, "features.rds"),
        greylist =  os.path.join(grey_path, "greylists.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        hic = os.path.join(annotation_path, "hic.rds"),
        regions = os.path.join(annotation_path,"gene_regions.rds"),
        results1 = os.path.join(
            diff_path, "{tgt1}", "{tgt1}_{comp1}-differential-signal.rds"
        ),
        results2 = os.path.join(
            diff_path, "{tgt2}", "{tgt2}_{comp2}-differential-signal.rds"
        ),
        script = os.path.join("workflow", "scripts", "pairwise_results.R"),
        yaml = "config/params.yml"
    output:
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise-results.rds"
        ),
        bed = expand(
            os.path.join(
                pairs_path, "{{tgt1}}_{{comp1}}-{{tgt2}}_{{comp2}}",
                "{{tgt1}}_{{comp1}}-{{tgt2}}_{{comp2}}-{f}.bed.gz"
            ),
            f = pw_dirs
        )
    params:
        pairwise_params = lambda wildcards: pairwise_param[wildcards.tgt1 + "_" + wildcards.comp1 + "-" + wildcards.tgt2 + "_" + wildcards.comp2] 
    threads: 2
    conda: "../envs/rmarkdown.yml"
    resources:
        runtime = "15m",
        mem_mb = 16000
    log: os.path.join(log_path, "prepare_pairwise_results", "{tgt1}_{comp1}-{tgt2}_{comp2}.log")
    script:
        "../scripts/pairwise_results.R"
    