localrules:
    generate_gene_query,
    get_protein_db,
    get_indexed_db,


rule get_indexed_db:
    output:
        "protein_annotation/index/Uniref50.lba.gz",
    params:
        "",
    log:
        "logs/lambda/get_indexed_db.log",
    conda:
        "../envs/base.yml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        wget -nv -q -O Uniref50.lba.gz http://ftp.imp.fu-berlin.de/pub/lambda/index/lambda3/gen_0/uniref50_20230713.lba.gz && \
        mv Uniref50.lba.gz {output} 2> {log}
        """


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
        num_matches=f'{config["lambda"]["num_matches"]}',
    log:
        "logs/lambda/blast_genes.log",
    conda:
        "../envs/lambda3.yml"
    shell:
        "lambda3 searchp -q {input.query} -i {input.indexed_db} -o {output} -n {params.num_matches} 2> {log}"
