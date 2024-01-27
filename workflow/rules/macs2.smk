def get_merged_bam_from_treat_and_target(wildcards):
    ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
    return expand(
        os.path.join(bam_path, "{file}.bam"), file = set(df[ind]['sample'])
    )

def get_merged_bai_from_treat_and_target(wildcards):
    ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
    return expand(
        os.path.join(bam_path, "{file}.bam.bai"), file = set(df[ind]['sample'])
    )

def get_input_bam_from_treat_and_target(wildcards):
    ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
    return expand(
        os.path.join(bam_path, "{file}.bam"),
        file = set(df[ind]['input'])
    )

def get_input_bai_from_treat_and_target(wildcards):
    ind = (df.treat == wildcards.treat) & (df.target == wildcards.target)
    return expand(
        os.path.join(bam_path, "{file}.bam.bai"),
        file = set(df[ind]['input'])
    )

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
    log: log_path + "/macs2_individual/{sample}.log"
    conda: "../envs/macs2.yml"
    params:
        outdir = os.path.join(macs2_path, "{sample}"),
        prefix = "{sample}",
        gsize = config['peaks']['macs2']['gsize'],
        fdr = config['peaks']['macs2']['fdr'],
        keep_duplicates = config['peaks']['macs2']['keep_duplicates']
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
            -n {params.prefix} \
            --bdg --SPMR \
            --outdir {params.outdir} 2> {log}
        cp {log} {output.log}
        """

rule macs2_qc:
    input:
        aln = lambda wildcards: expand(
            os.path.join(bam_path, "{sample}.{suffix}"),
            sample = set(df[df.target == wildcards.target]['sample']),
            suffix = ['bam', 'bam.bai']
        ),
        blacklist = blacklist,
        config = "config/config.yml",
        extrachips = rules.update_extrachips.output,
        indiv_macs2 = lambda wildcards: expand(
            os.path.join(macs2_path, "{sample}", "{sample}_{suffix}"),
            sample = set(df[df.target == wildcards.target]['sample']),
            suffix = ['callpeak.log', 'peaks.narrowPeak']
        ),
        samples = config['samples']['file'],
        seqinfo = os.path.join(annotation_path, "seqinfo.rds"),
        r = "workflow/scripts/macs2_qc.R"
    output:
        cors = os.path.join(
            macs2_path, "{target}", "{target}_cross_correlations.tsv"
        ),
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv")
    conda: "../envs/rmarkdown.yml"
    threads: lambda wildcards: len(df[df['target'] == wildcards.target])
    resources:
        mem_mb = 8192
    log: log_path + "/macs2_qc/{target}_macs2_qc.log"
    shell:
        """
        ## Run the QC script
        Rscript --vanilla \
            {input.r} \
            {wildcards.target} \
            {threads} &>> {log}
        """

rule macs2_merged:
    input:
        bam = get_merged_bam_from_treat_and_target,
        bai = get_merged_bai_from_treat_and_target,
        control = get_input_bam_from_treat_and_target,
        control_bai = get_input_bai_from_treat_and_target,
        qc = os.path.join(macs2_path, "{target}", "{target}_qc_samples.tsv")
    output:
        narrow_peaks = os.path.join(
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
    log: log_path + "/macs2_merged/{target}/{treat}_merged.log"
    conda: "../envs/macs2.yml"
    params:
        bamdir = bam_path,
        outdir = os.path.join(macs2_path, "{target}"),
        prefix = "{target}_{treat}_merged",
        gsize = config['peaks']['macs2']['gsize'],
        fdr = config['peaks']['macs2']['fdr'],
        keep_duplicates = config['peaks']['macs2']['keep_duplicates']
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
    log:
        "workflow/logs/macs2_bdgcmp/{target}_{treat}_bdgcmp.log"
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
        
