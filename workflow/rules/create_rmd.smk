def get_bw_type(wildcards):
    heat_type = profile_heatmap_param[wildcards.target]['bw_type']
    if heat_type == 'FE':
        bw_type = ['FE']
    else:
        bw_type = ['treat_pileup']
    path = expand(
        os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_{bw_type}.bw"
        ),
        bw_type = bw_type,
        target = [wildcards.target], treat = [wildcards.ref, wildcards.treat]
    )
    return(path)



rule create_site_yaml:
    input:
        arg_checks = rules.check_args.output,
        packages = rules.check_r_packages.output,
        samples = config['samples']['file'],
        script = os.path.join("workflow", "scripts", "create_site_yaml.R"),
        yml = "config/rmarkdown.yml",
    output:
        yml = os.path.join(rmd_path, "_site.yml")
    log: os.path.join(log_path, "create_rmd", "site_yaml.log")
    threads: 1
    localrule: True
    resources:
        mem_mb = 1024,
        runtime = "5m",
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/create_site_yaml.R"

rule create_setup_chunk:
    input:
        packages = rules.check_r_packages.output,
        script = os.path.join("workflow", "scripts", "create_setup_chunk.R"),
        yml = "config/rmarkdown.yml",
    output:
        rmd = "analysis/setup_chunk.Rmd"
    log: os.path.join(log_path, "create_rmd", "setup_chunk.log")
    threads: 1
    localrule: True
    resources:
        mem_mb = 1024,
        runtime = "5m",
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/create_setup_chunk.R"

rule create_index_rmd:
    input:
        packages = rules.check_r_packages.output,
        rmd = os.path.join("workflow", "modules", "index.Rmd"),
    output:
        os.path.join(rmd_path, "index.Rmd")
    threads: 1
    localrule: True
    resources:
        mem_mb = 512,
        runtime = "2m",
    shell:
        """
        cat {input.rmd} > {output}
        """

rule create_annotations_rmd:
    input:
        arg_checks = rules.check_args.output,
        blacklist = blacklist,
        chrom_sizes = chrom_sizes,
        features = rules.prep_features.output.rds, 
        gene_regions = rules.create_genome_annotations.output.regions, 
        greylist = rules.combine_greylists.output.rds,
        gsea_dir = os.path.join(annotation_path, "gsea_dir.rds"),
        gsea_sig = os.path.join(annotation_path, "gsea_sig.rds"),
        gtf = rules.create_genome_annotations.output.gtf,
        hic = rules.prep_hic.output.hic, 
        module = os.path.join(
            "workflow", "modules", "annotation_description.Rmd"
        ),
        motifs = rules.prep_motifs.output.motifs,
        motif_uri = rules.prep_motifs.output.motif_uri,    
        packages = rules.check_r_packages.output,   
        rna = os.path.join(annotation_path, "rna.rds"),    
        script = os.path.join(
            "workflow", "scripts", "create_annotations_rmd.R"
        ),
        seqinfo = rules.create_genome_annotations.output.seqinfo, 
        trans_models = os.path.join(annotation_path, "trans_models.rds"),
        tss = os.path.join(annotation_path, "tss.rds"),
    output:
        rmd = os.path.join(rmd_path, "annotation_description.Rmd"),
    params:
        colours = os.path.join(annotation_path, "colours.rds"),
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "create_rmd", "annotation_description.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_annotations_rmd.R"

rule create_signal_summary_rmd:
    input:
        arg_checks = rules.check_args.output,
        bw = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}",
                "{{target}}_{treat}_merged_treat_pileup.bw"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        cors = os.path.join(
            macs2_path, "{target}", "{target}_cross_correlations.tsv"
        ),
        module = "workflow/modules/signal_summary.Rmd",
        packages = rules.check_r_packages.output,
        peak_files = expand(
             os.path.join(
                peak_path, "{{target}}", "{{target}}_{f}"
            ),
            f = [
                'motif_position.tsv.gz', 'motif_enrichment.tsv.gz',
                'regions_localz.rds', 'consensus_peaks.bed.gz'
                ]
        ),
        script = os.path.join("workflow", "scripts", "create_signal_summary.R"),
    output:
        rmd = os.path.join(rmd_path, "{target}_signal_summary.Rmd")
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "create_rmd", "create_{target}_signal_summary.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_signal_summary.R"


rule create_nfr_rmd:
    input:
        arg_checks = rules.check_args.output,
        bigwig = lambda wildcards: expand(
            os.path.join(
                macs2_path, "{{target}}", 
                "{{target}}_{treat}_merged_treat_pileup.bw"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        consensus_peaks = expand(
            os.path.join(peak_path, "{t}", "{t}_consensus_peaks.rds"),
            t = targets
        ),
        files = expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_consensus_nfr.{suffix}"
            ),
            suffix = ['rds', 'bed.gz']
        ),
        localz = expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_nfr_{f}_localz.rds"
            ),
            f = ['regions', 'targets']
        ),
        module = os.path.join("workflow", "modules", "nfr.Rmd"),
        motif_results = expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_motif_{f}.tsv.gz"
            ),
            f = ['enrichment', 'position']
        ),
        nfr = lambda wildcards: expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_{treat}.nfr.bed.gz"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        packages = rules.check_r_packages.output,
        peaks = lambda wildcards: expand(
            os.path.join(
                nfr_path, "{{target}}", "{{target}}_{treat}_merged_peaks.bed"
            ),
            treat = set(df[df.target == wildcards.target]['treat'])
        ),
        script = os.path.join("workflow", "scripts", "create_nfr_rmd.R"),
    output:
        rmd = os.path.join(rmd_path, "{target}_nfr.Rmd")
    conda: "../envs/rmarkdown.yml"
    threads: 1
    localrule: True
    log: os.path.join(log_path, "create_rmd", "create_{target}_nfr.log")
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_nfr_rmd.R"

rule create_differential_signal_rmd:
    input:
        annotations = ANNOTATION_RDS,
        arg_checks = rules.check_args.output,
        bigwig = get_bw_type,
        counts = os.path.join(diff_path, "{target}", "{target}_counts.rds"),
        ihw = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-ihw.rds"
        ),
        localz = os.path.join(
            diff_path, "{target}", "{target}_{ref}_{treat}-regions_localz.rds"
        ),
        module = os.path.join("workflow", "modules", "differential_signal.Rmd"),
        motif_enrichment = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}_motif_enrichment.tsv.gz"
        ),
        motif_position = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}_motif_position.tsv.gz"
        ),
        nfr = NFR_RDS,
        packages = rules.check_r_packages.output,
        r = os.path.join("workflow", "scripts", "create_differential_rmd.R"),
        results = os.path.join(
            diff_path, "{target}", 
            "{target}_{ref}_{treat}-differential-signal.rds"
        )
    output:
        rmd = os.path.join(
            rmd_path, "{target}_{ref}_{treat}_differential_signal.Rmd"
        )
    conda: "../envs/rmarkdown.yml"
    localrule: True
    log: os.path.join(log_path, "create_rmd", "{target}_{ref}_{treat}_differential_signal.log")
    threads: 1
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_differential_rmd.R"


rule create_pairwise_comparisons_rmd:
    input:
        arg_checks = rules.check_args.output,
        dsa1 = os.path.join(
            diff_path, "{tgt1}", "{tgt1}_{comp1}-differential-signal.rds"
        ),
        dsa2 = os.path.join(
            diff_path, "{tgt2}", "{tgt2}_{comp2}-differential-signal.rds"
        ),      
        localz = rules.merge_localz_pairwise.output.rds,
        module = os.path.join(
            "workflow", "modules", "pairwise_comparison.Rmd"
        ),
        motif_enrich = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-motif_enrichment.tsv.gz"
        ),
        motif_position = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-motif_position.tsv.gz"
        ),   
        packages = rules.check_r_packages.output,
        results = os.path.join(
            pairs_path, "{tgt1}_{comp1}-{tgt2}_{comp2}", 
            "{tgt1}_{comp1}-{tgt2}_{comp2}-pairwise-results.rds"
        ),
        script = os.path.join("workflow", "scripts", "create_pairwise_rmd.R")
    output:
        rmd = os.path.join(
            rmd_path, "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise_comparison.Rmd"
        )
    conda: "../envs/rmarkdown.yml"
    localrule: True
    log: os.path.join(log_path, "create_rmd", "{tgt1}_{comp1}-{tgt2}_{comp2}_pairwise.log")
    threads: 1
    resources:
        mem_mb = 1024,
        runtime = "5m",
    script:
        "../scripts/create_pairwise_rmd.R"