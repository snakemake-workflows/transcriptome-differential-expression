localrules:
    get_genome,
    extract_annotation,
    extract_genome,


rule get_genome:
    output:
        # generic name:
        temp("references/ncbi_dataset.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_genome.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include gff3,genome &> {log}
        """


rule extract_genome:
    input:
        lambda wildcards: get_reference_files(config).get("genome","references/ncbi_dataset.zip")
    output:
        "references/genomic.fa",
    group:
        "reference"
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/extract_genome.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        [ -f "{input}" ] && mv {input} {output} || unzip -p {input} ncbi_dataset/data/{params.accession}/*.fna > {output} 2>> {log}
        """


rule extract_annotation:
    input:
        lambda wildcards: get_reference_files(config).get("annotation","references/ncbi_dataset.zip")
    output:
        "references/genomic.gff",
    group:
        "reference"
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        [ -f "{input}" ] && mv {input} {output} || unzip -p {input} ncbi_dataset/data/{params.accession}/*.gff > references/genomic.gff 2>> {log};
        """
