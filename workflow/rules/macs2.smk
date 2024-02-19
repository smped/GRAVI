
rule macs2_individual:
    input:
        bam = os.path.join(bam_path, "{sample}.bam"),
        bai = os.path.join(bam_path, "{sample}.bam.bai"),
        control = lambda wildcards: expand(
            os.path.join(bam_path, "{ctrl}.{suffix}"),
            ctrl = set(df[df['sample'] == wildcards.sample]['input']),
            suffix = ['bam', 'bam.bai']
        ),
    output:
        narrow_peaks = os.path.join(
            macs2_path, "{sample}", "{sample}_peaks.narrowPeak"
        ),
        bedgraph = temp(
            os.path.join(macs2_path, "{sample}", "{sample}_treat_pileup.bdg")
        ),
        log = os.path.join(macs2_path, "{sample}", "{sample}_callpeak.log"),
        summits = os.path.join(macs2_path, "{sample}", "{sample}_summits.bed"),
        other = temp(
            expand(
                os.path.join(macs2_path, "{{sample}}", "{{sample}}{suffix}"),
                suffix = ['_model.r', '_peaks.xls', '_control_lambda.bdg']
            )
        )
    log: os.path.join(log_path, "macs2_individual", "{sample}.log")
    conda: "../envs/macs2.yml"
    params:
        outdir = os.path.join(macs2_path, "{sample}"),
        gsize = lambda wildcards: macs2_param[
            list(df[df['sample'] == wildcards.sample]['target'])[0]
        ]['gsize'],
        fdr = lambda wildcards: macs2_param[
            list(df[df['sample'] == wildcards.sample]['target'])[0]
        ]['fdr'],
        keep_duplicates = lambda wildcards: macs2_param[
            list(df[df['sample'] == wildcards.sample]['target'])[0]
        ]['keep_duplicates']
    threads: 1
    resources:
        mem_mb = 8192
    shell:
        """
        macs2 callpeak \
            -t {input.bam}\
            -c {input.control[0]} \
            -f BAM \
            -g {params.gsize} \
            --keep-dup {params.keep_duplicates} \
            -q {params.fdr} \
            -n {wildcards.sample} \
            --bdg --SPMR \
            --outdir {params.outdir} 2> {log}
        cp {log} {output.log}
        """

rule macs2_qc:
    input:
        bam = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.bam"),
            sample = set(df[df.target == wildcards.target]['sample']),
        ),
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        checks = ALL_CHECKS,
        greylist = lambda wildcards: expand(
            os.path.join(annotation_path, "{f}_greylist.bed.gz"),
            f = set(df[df.target == wildcards.target]['input'])
        ),
        input_bam = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.bam"),
            sample = set(df[df.target == wildcards.target]['input']),
        ),
        logs = lambda wildcards: expand(
            os.path.join(macs2_path, "{sample}", "{sample}_callpeak.log"),
            sample = set(df[df.target == wildcards.target]['sample']),
        ),
        peaks = lambda wildcards: expand(
            os.path.join(macs2_path, "{sample}", "{sample}_peaks.narrowPeak"),
            sample = set(df[df.target == wildcards.target]['sample']),
        ),
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
        gene_regions = os.path.join(annotation_path, "gene_regions.rds"),
    output:
        cors = os.path.join(
            macs2_path, "{target}", "{target}_cross_correlations.tsv"
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv"),
    params:
        outlier_threshold = lambda wildcards: macs2_qc_param[wildcards.target]['outlier_threshold'],
        allow_zero = lambda wildcards: macs2_qc_param[wildcards.target]['allow_zero'],
    conda: "../envs/rmarkdown.yml"
    threads: lambda wildcards: len(df[df['target'] == wildcards.target])
    resources:
        mem_mb = 16384,
    retries: 1
    log: os.path.join(log_path, "macs2_qc", "{target}_macs2_qc.log")
    script:
        "../scripts/macs2_qc.R"

rule macs2_merged:
    input:
        bam = lambda wildcards: expand(
            os.path.join(bam_path, "{f}.bam"),
            f = set(
                df[(df.treat == wildcards.treat) & (df.target == wildcards.target)]['sample']
            )
        ),
        bai = lambda wildcards: expand(
            os.path.join(bam_path, "{f}.bam.bai"),
            f = set(
                df[(df.treat == wildcards.treat) & (df.target == wildcards.target)]['sample']
            )
        ),
        input = lambda wildcards: expand(
            os.path.join(bam_path, "{f}.{suffix}"),
            f = set(
                df[(df.treat == wildcards.treat) & (df.target == wildcards.target)]['input']
            ),
            suffix = ['bam', 'bam.bai']
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv")
    output:
        peaks = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_peaks.narrowPeak"
        ),
        summits = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_summits.bed"
        ),
        bedgraph = temp(
            expand(
                os.path.join(
                    macs2_path, "{{target}}", "{{target}}_{{treat}}_merged_{type}.bdg"
                ),
                type = ['treat_pileup', 'control_lambda']
            )
        ),
        log = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_callpeak.log"
        )
    log: os.path.join(log_path, "macs2_merged", "{target}_{treat}_merged.log")
    conda: "../envs/macs2.yml"
    params:
        bamdir = bam_path,
        outdir = os.path.join(macs2_path, "{target}"),
        prefix = "{target}_{treat}_merged",
        gsize = lambda wildcards: macs2_param[wildcards.target]['gsize'],
        fdr = lambda wildcards: macs2_param[wildcards.target]['fdr'],
        keep_duplicates = lambda wildcards: macs2_param[wildcards.target]['keep_duplicates']
    threads: 1
    resources:
        mem_mb = 8192	
    shell:
        """
        QC_PASS=$(egrep 'pass' {input.qc} | egrep {wildcards.treat} | cut -f1)
        SAMPLES=$(for f in $QC_PASS; do echo {params.bamdir}/$f.bam; done)

        ## Get the input column
        I=$(head -n1 {input.qc} | sed -r 's/\\t/\\n/g' | egrep -n '[Ii]nput' | sed -r 's/([0-9]+):[Ii]nput/\\1/g')
        INPUT_PASS=$(egrep 'pass$' {input.qc} | egrep "{wildcards.treat}\s" | cut -f$I | uniq)
        INPUT=$(for f in $INPUT_PASS; do echo {params.bamdir}/$f.bam; done)

        macs2 callpeak \
            -t $SAMPLES\
            -c $INPUT \
            -f BAM \
            -g {params.gsize} \
            --keep-dup {params.keep_duplicates} \
            -q {params.fdr} \
            -n {params.prefix} \
            --bdg --SPMR \
            --outdir {params.outdir} 2> {log}
        cp {log} {output.log}
        """

rule macs2_bdgcmp:
    input:
        bdg = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_treat_pileup.bdg"
        ),
        ctrl = os.path.join(
            macs2_path, "{target}", "{target}_{treat}_merged_control_lambda.bdg"
        ),
    output:
        temp(
            os.path.join(macs2_path, "{target}/{target}_{treat}_merged_FE.bdg"),
        )
    log: os.path.join(log_path, "macs2_bdgcmp", "{target}", "{target}_{treat}_bdgcmp.log")
    conda: "../envs/macs2.yml"
    threads: 1
    resources:
        mem_mb = 16384,
        runtime = "3h"
    shell:
        """
        macs2 bdgcmp \
            -t {input.bdg} \
            -c {input.ctrl} \
            -m FE \
            -o {output} 2> {log}
        """		
        
rule filter_merged_peaks:
    input:
        blacklist = os.path.join(annotation_path, "blacklist.rds"),
        greylist = lambda wildcards: expand(
            os.path.join(annotation_path, "{f}_greylist.bed.gz"),
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
            macs2_path, "{target}", "{target}_{treat}_filtered_peaks.narrowPeak"
        )
    params:
        min_prop = lambda wildcards: macs2_qc_param[wildcards.target]['min_prop_reps']
    conda: "../envs/rmarkdown.yml"
    threads: 1
    log: os.path.join(log_path, "filter_merged_peaks", "{target}_{treat}.log")
    resources:
        mem_mb = 4096,
        runtime = "15m"
    script:
        "../scripts/filter_merged_peaks.R"
