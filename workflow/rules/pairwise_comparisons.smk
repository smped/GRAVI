rule prepare_pairwise_results:
    input:
        arg_checks = rules.check_args.output,
        packages = rules.check_r_packages.output,    
        blacklist = rules.prep_blacklist.output.blacklist, 
        features = rules.prep_features.output.rds, 
        greylist =  rules.combine_greylists.output.rds,
        gtf = rules.create_genome_annotations.output.gtf,
        hic = rules.prep_hic.output.hic, 
        regions = rules.create_genome_annotations.output.regions, 
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

rule motif_analysis_pairwise:
    input:
        arg_checks = rules.check_args.output,
        rds = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise-results.rds"
        ),
        motifs = rules.prep_motifs.output.motifs,
        packages = rules.check_r_packages.output,
        script = os.path.join(
            "workflow", "scripts", "pairwise_motif_analysis.R"
        ),
        seqinfo = rules.create_genome_annotations.output.seqinfo, 
    output:
        enrich_tsv = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-motif_enrichment.tsv.gz"
        ),
        position_tsv = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-motif_position.tsv.gz"
        ),        
    params:
        motif_params = motif_param['pairwise']
    retries: 2
    threads: lambda wildcards, attempt: attempt * 8
    conda: "../envs/rmarkdown.yml"
    resources:
        disk_mb = lambda wildcards, attempt: attempt * 8000,
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        runtime = lambda wildcards, attempt: attempt * 30,
    log: os.path.join(log_path, "pairwise_motif_analysis", "{tgt1}_{comp1}-{tgt2}_{comp2}.log")
    script:
        "../scripts/pairwise_motif_analysis.R"    

