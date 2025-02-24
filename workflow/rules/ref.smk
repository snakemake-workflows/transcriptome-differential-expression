localrules:
    get_references,
    get_annotation,
    get_genome,


rule get_references:
    output:
        # we need two different output, to ensure simultaneous access
        # by the two downstream rules:
        a=temp("references/ncbi_dataset_a.zip"),
        b=temp("references/ncbi_dataset_b.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_references.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include gff3,genome &> {log} && mv ncbi_dataset.zip {output.a};
        cp {output.a} {output.b}
        """


rule get_genome:
    input:
        lambda wildcards: get_reference_files(config).get("genome"),
    output:
        temp("references/genomic.fa"),
    priority: 10
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_genome.log",
    conda:
        "../envs/reference.yml"
    script:
        "../scripts/extract_refs.py"


rule get_annotation:
    input:
        lambda wildcards: get_reference_files(config).get("annotation"),
    output:
        temp("references/genomic.gff"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/reference.yml"
    script:
        "../scripts/extract_refs.py"
