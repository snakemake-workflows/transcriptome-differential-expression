rule genome_to_transcriptome:
    input:
        genome=config["genome"],
        annotation=config["annotation"],
    output:
        transcriptome="transcriptome/transcriptome.fa",
    log:
        "logs/gffread.log",
    conda:
        "../envs/env.yml"
    shell:
        """
    gffread -t {resources.cpus_per_task} -w {output} -g {input.genome} {input.annotation} &> {log}    
    """


rule filter_reads:
    input:
        fastq=lambda wildcards: get_mapped_reads_input(
            samples["sample"][wildcards.sample]
        ),
    output:
        "filter/{sample}_filtered.fq",
    message:
        f"Filtering with read length >= {config['min_length']}."
    log:
        "logs/filter_reads/{sample}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/read_filter.py"
