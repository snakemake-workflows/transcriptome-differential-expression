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
    shell:
        """
        curl -o data.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.accession}/download?include_annotation_type=GENOME_FASTA &> {log};
        unzip -p data.zip ncbi_dataset/data/{params.accession}/*.fna > references/genomic.fa 2> {log};
        rm data.zip &> {log}
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
    shell:
        """
        curl -o data.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.accession}/download?include_annotation_type=GENOME_GFF &> {log};
        unzip -p data.zip ncbi_dataset/data/{params.accession}/*.gff > references/genomic.gff 2> {log};
        rm data.zip &> {log}
        """
