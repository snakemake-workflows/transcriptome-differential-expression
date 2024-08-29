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
        "../envs/curl.yml"
    shell:
        """
        curl -s -o data_genome.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.accession}/download?include_annotation_type=GENOME_FASTA &> {log};
        unzip -p data_genome.zip ncbi_dataset/data/{params.accession}/*.fna > references/genomic.fa 2> {log};
        rm data_genome.zip &> {log}
        """


rule get_annotation:
    output:
        "references/genomic.gff",
    params:
        accession=config["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/curl.yml"
    shell:
        """
        curl -s -o data_annotation.zip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.accession}/download?include_annotation_type=GENOME_GFF &> {log};
        unzip -p data_annotation.zip ncbi_dataset/data/{params.accession}/*.gff > references/genomic.gff 2> {log};
        rm data_annotation.zip &> {log}
        """
