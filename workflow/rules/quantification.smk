localrules:
    merge_counts,


rule count_reads:
    input:
        bam="sorted_alignments/{sample}.bam",
        trs="transcriptome/transcriptome.fa",
    output:
        tsv=temp("counts/{sample}_salmon/quant.sf"),
    params:
        outdir=temp(lambda wildcards: f"counts/{wildcards.sample}_salmon"),
        libtype=config["salmon_libtype"],
    log:
        "logs/salmon/{sample}.log",
    conda:
        "../envs/salmon.yml"
    shell:
        """
        salmon --no-version-check quant --ont -p {resources.cpus_per_task} \
        -t {input.trs} -l {params.libtype} -a {input.bam} -o {params.outdir} 2> {log}
        """


rule merge_counts:
    input:
        count_tsvs=expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"]),
    output:
        temp("merged/all_counts.tsv"),
    log:
        "logs/merge_count.log",
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/merge_count_tsvs.py"
