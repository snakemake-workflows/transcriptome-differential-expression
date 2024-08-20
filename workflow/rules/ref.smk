localrules:
    get_genome,
    get_annotation,


rule get_genome:
    output:
        genome="references/genomic.fa",
    params:
        accession=config["accession"],
    log:
        "logs/refs/get_genome.log",
    conda:
        "../envs/env.yml"
    script:
        """
        curl -o data.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.accession}/download?include_annotation_type=GENOME_FASTA;
        mkdir -p references;
        unzip -p data.zip ncbi_dataset/data/{params.accession}/*.fna > {output.annotation};
        rm data.zip;
        """


rule get_annotation:
    output:
        "references/genomic.gff",
    params:
        accession=config["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/env.yml"
    script:
        """
        mkdir -p references
        curl -o data.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.accession}/download?include_annotation_type=GENOME_GFF
        && unzip -p data.zip ncbi_dataset/data/{params.accession}/*.gff > references/genomic.gff &&
        && rm data.zip
        """
