localrules:
    generate_gene_query,
    get_protein_db,


rule get_protein_db:
    output:
        "protein_annotation/index/Uniref50.fa.gz",
    params:
        "",
    log:
        "logs/lambda/get_protein_db.log",
    conda:
        "../envs/base.yml"
    shell:
        """
        wget -nv -O {output}  ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
        &> {log}
        """


rule index_protein_db:
    input:
        "protein_annotation/index/Uniref50.fa.gz",
    output:
        "protein_annotation/index/Uniref50.lba.gz",
    params:
        "",
    log:
        "logs/lambda/index_protein_db.log",
    conda:
        "../envs/lambda3.yml"
    shell:
        "lambda3 mkindexp -d {input} 2> {log}"


rule generate_gene_query:
    input:
        sorted_lfc_counts="de_analysis/sorted_normalized_counts.csv",
        transcriptome="transcriptome/transcriptome.fa",
    output:
        "protein_annotation/de_genes.fa",
    params:
        "",
    log:
        "logs/lambda/generate_gene_query.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/get_de_genes.py"


rule blast_genes:
    input:
        indexed_db="protein_annotation/index/Uniref50.lba.gz",
        query="protein_annotation/de_genes.fa",
    output:
        "protein_annotation/blast_results.m8",
    params:
        "",
    log:
        "logs/lambda/blast_genes.log",
    conda:
        "../envs/lambda3.yml"
    shell:
        "lambda3 searchp -q {input.query} -i {input.indexed_db} -o {output} 2> {log}"
