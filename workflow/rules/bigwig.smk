rule sort_bedgraph:
    input: 
        os.path.join(macs2_path, "{path}", "{f}.bdg"),
    output:
        temp(os.path.join(macs2_path, "{path}", "{f}.sorted.bdg"))
    log: os.path.join(log_path, "sort_bedgraph", "{path}", "{f}.log")
    threads: 2
    params:
        pattern = "$'^chr[0-9XY]+\\t'" # alt: "$'^[0-9XY]+\\t'"
    resources:
        runtime = "1h",
        mem_mb = lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8,
        disk_mb = lambda wildcards, input, attempt: (input.size//1000000) * attempt * 4,
    shell:
        """
        ## Sort the file
        echo -e "Started sorting at $(date)" >> {log}
        sort \
          -k1,1 -k2,2n \
          -S {resources.mem_mb}M \
          --parallel {threads} \
          {input} | \
          egrep {params.pattern} > {output}
        echo -e "Finished sorting at $(date)" >> {log}
        """	


rule bedgraph_to_bigwig:
    input:
        bedgraph = os.path.join(macs2_path, "{path}", "{f}.sorted.bdg"),
        chrom_sizes = chrom_sizes
    output:
        bigwig = os.path.join(macs2_path, "{path}", "{f}.bw")
    conda: "../envs/bedgraph_to_bigwig.yml"
    log: os.path.join(log_path, "bedgraph_to_bigwig", "{path}", "{f}.log")
    # This gives ~1 thread per 500Mb & should track with mem_mb resources
    threads: 8
    resources:
        runtime = "2h",
        mem_mb = lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8,
        disk_mb = lambda wildcards, input, attempt: (input.size//1000000) * attempt * 4,
    shell:
        """
        echo -e "Started conversion at $(date)" >> {log}
        bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bigwig}
        echo -e "Finished conversion at $(date)" >> {log}
        """

