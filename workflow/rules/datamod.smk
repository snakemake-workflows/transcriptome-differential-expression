localrules:
    genome_to_transcriptome,


rule genome_to_transcriptome:
    input:
        genome="references/genomic.fa",
        annotation="references/genomic.gff",
    output:
        transcriptome=temp("transcriptome/transcriptome.fa"),
    log:
        "logs/gffread.log",
    conda:
        "../envs/gffread.yml"
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
        temp("filter/{sample}_filtered.fq"),
    message:
        f"Filtering with read length >= {config['min_length']}."
    log:
        "logs/filter_reads/{sample}.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/read_filter.py"
