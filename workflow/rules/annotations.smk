rule create_genome_annotations:
    input:
        bam = expand(os.path.join(bam_path, "{bam}.bam"), bam = samples),
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        script = os.path.join("workflow", "scripts", "create_genome_annotations.R"),
        yaml = os.path.join("config", "params.yml"),
    output:
        gtf_exon = os.path.join(annotation_path, "gtf_exon.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        gtf_transcript = os.path.join(annotation_path, "gtf_transcript.rds"),
        regions = os.path.join(annotation_path, "gene_regions.rds"),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
        trans_models = os.path.join(annotation_path, "trans_models.rds"),
        tss = os.path.join(annotation_path, "tss.rds"),
    conda: "../envs/rmarkdown.yml"
    params:
        colours = os.path.join(annotation_path, "colours.rds"),
    threads: 2
    retries: 1
    resources:
        mem_mb = 16384,
        run_time = "30m"
    log: os.path.join(log_path, "annotations", "genome_annotations.log")
    script:
        "../scripts/create_genome_annotations.R"

## These might need to be renamed as they're the N-rich regions rather than
## blacklisted. Very useful for motif analysis
rule make_exclude_ranges:
    input:
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        script = os.path.join(
            "workflow", "scripts", "make_exclude_ranges.R"
        ),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds")
    output:
        rds = os.path.join(annotation_path, "exclude_ranges.rds")
    threads: 1
    localrule: True
    retries: 1    
    resources:
        runtime = "10m",
        mem_mb = 4096,
    log: os.path.join(log_path, "annotations", "make_exclude_ranges.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/make_exclude_ranges.R"

rule prep_features:
    input:
        gene_regions = os.path.join(annotation_path, "gene_regions.rds"),
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),    
        script = os.path.join("workflow", "scripts", "prep_features.R"),
    output:
        rds = os.path.join(annotation_path, "features.rds"),
    threads: 1
    retries: 1    
    resources:
        mem_mb = 8192,
        run_time = "20m"
    log: os.path.join(log_path, "annotations", "features.log")
    conda: "../envs/rmarkdown.yml"
    script:
        "../scripts/prep_features.R"

rule prep_blacklist:
    input:
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
        script = os.path.join("workflow", "scripts", "prep_blacklist.R"),
    output:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1    
    resources:
        mem_mb = 8192,
        run_time = "10m"
    log: os.path.join(log_path, "annotations", "blacklist.log")
    script:
        "../scripts/prep_blacklist.R"


rule prep_hic:
    input:
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
        script = os.path.join("workflow", "scripts", "prep_hic.R"),
        yaml = os.path.join("config", "params.yml"),
    output:
        hic = os.path.join(annotation_path, "hic.rds"),
    conda: "../envs/rmarkdown.yml"
    threads: 2
    retries: 1    
    resources:
        mem_mb = 16000,
        run_time = "10m"
    log: os.path.join(log_path, "annotations", "hic.log")
    script:
        "../scripts/prep_hic.R"

rule prep_motifs:
    input:
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        script = os.path.join("workflow", "scripts", "prep_motifs.R"),
        yaml = os.path.join("config", "params.yml"),
    output:
        motifs = os.path.join(annotation_path, "motif_list.rds"),
        motif_uri = os.path.join(annotation_path, "motif_uri.rds"),
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1    
    resources:
        mem_mb = 8192,
        run_time = "10m"
    log: os.path.join(log_path, "annotations", "prep_motifs.log")
    script:
        "../scripts/prep_motifs.R"

rule prep_msigdb:
    input:
        here = os.path.join(check_path, "here.chk"),
        packages = os.path.join(check_path, "r-packages.chk"),
        script = os.path.join("workflow", "scripts", "prep_msigdb.R"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        yaml = os.path.join("config", "params.yml"),
    output:
        msigdb = os.path.join(annotation_path, "msigdb.rds"),
    conda: "../envs/rmarkdown.yml"
    threads: 1
    retries: 1    
    resources:
        mem_mb = 8192,
        run_time = "10m"
    log: os.path.join(log_path, "annotations", "msigdb.log")
    script:
        "../scripts/prep_msigdb.R"

rule prep_rna:
    input:
        here = os.path.join(check_path, "here.chk"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        msigdb = os.path.join(annotation_path, "msigdb.rds"),
        packages = os.path.join(check_path, "r-packages.chk"),
        script = os.path.join("workflow", "scripts", "prep_rna.R"),
        yaml = os.path.join("config", "params.yml"),
    output:
        gsea_dir = os.path.join(annotation_path, "gsea_dir.rds"),
        gsea_sig = os.path.join(annotation_path, "gsea_sig.rds"),
        rna = os.path.join(annotation_path, "rna.rds"),
    params:
        gene_col = ["gene_id", "Geneid", "geneid", "ensembl_gene_id", "ensembl_id"],
        expr_col = ["AveExpr", "logCPM", "baseMean"],
        lfc_col = ["logFC", "logfc", "lfc", "log2FoldChange"],
        p_col = ["PValue", "PVal", "P", "p", "p_value", "p_val", "P.Value", "pvalue"],
        padj_col = ["fdr", "FDR", "adjP", "adj_p", "adj.P.Value", "padj"],
        nperm_gsea = 1e5,
    conda: "../envs/rmarkdown.yml"
    threads: 4
    retries: 1    
    resources:
        mem_mb = 16384,
        run_time = "30m"
    log: os.path.join(log_path, "annotations", "rna.log")
    script:
        "../scripts/prep_rna.R"

rule make_chrom_sizes:
    input: expand(os.path.join(bam_path, "{sample}.bam"), sample = [samples[0]])
    output: chrom_sizes
    conda: "../envs/samtools.yml"
    threads: 1
    retries: 1    
    resources:
        run_time = "5m"
    shell:
        """
        samtools view -H {input} | \
          egrep '^@SQ' | \
          cut -f2,3 | \
          sed -r 's/^SN:(.+)\\tLN:(.+)$/\\1\\t\\2/g' |\
          egrep '^[c0-9]' |\
          egrep -v 'M' |\
          sort -V > {output}
        """