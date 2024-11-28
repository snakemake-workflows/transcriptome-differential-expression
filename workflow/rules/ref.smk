localrules:
    get_references,
    get_annotation,
    get_genome,


rule get_references:
    output:
        # generic name:
        temp("references/ncbi_dataset.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_references.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include gff3,genome &> {log} && mv ncbi_dataset.zip {output}
        """


rule get_genome:
    input:
        lambda wildcards: get_reference_files(config).get(
            "genome", "references/ncbi_dataset.zip"
        ),
    output:
        "references/genomic.fa",
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_genome.log",
    conda:
        "../envs/reference.yml"
    shell:
        # checks if local genome is available (see commons.smk), if it is, moves it to output path. If not, extract genome from ncbi_dataset.zip
        """
        [ -f "{input}" ] && cp "{input}" "{output}" 2> "{log}"  || unzip -p {input} ncbi_dataset/data/{params.accession}/*.fna > {output} 2>> {log}
        """


rule get_annotation:
    input:
        lambda wildcards: get_reference_files(config).get(
            "annotation", "references/ncbi_dataset.zip"
        ),
    output:
        "references/genomic.gff",
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/reference.yml"
    shell:
        # checks if local annotation is available (see commons.smk), if it is, moves it to output path. If not, extracts annotation from ncbi_dataset.zip
        """
        [ -f "{input}" ] && cp "{input}" "{output}" 2> "{log}" || unzip -p {input} ncbi_dataset/data/{params.accession}/*.gff > {output} 2>> {log};
        """
