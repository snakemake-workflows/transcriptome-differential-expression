localrules:
    samples_tab,


rule samples_tab:
    input:
        sorted_bams=expand(
            "sorted_alignments/{sample}_sorted.bam", sample=samples["sample"]
        ),
    output:
        "ESPRESSO/sorted_bam_samples.tsv",
    params:
        samples=samples["sample"],
    log:
        "logs/ESPRESSO/samples_tab.log",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/sorted_bam_tsv.py"


rule gff_to_gtf:
    input:
        "references/standardized_genomic.gff",
    output:
        "references/standardized_genomic.gtf",
    log:
        "logs/gffread/gff_to_gtf.log",
    conda:
        "../envs/gffread.yml"
    script:
        """
        gffread -t {resources.cpus_per_task} -E {input} -T -o {output} &> {log}    
        """

rule ESPRESSO_S:
    input:
        bams="ESPRESSO/sorted_bam_samples.tsv",
        transcriptome="transcriptome/transcriptome.fa",
        annotation="references/standardized_genomic.gtf",
    output:
        ""
    log:
        "logs/ESPRESSO_S.log"
    conda:
        "../envs/ESPRESSO.yml"
    script:
        ""


rule ESPRESSO_C:
    input:
    output:
    log:
        "logs/ESPRESSO_C.log"
    conda:
        "../envs/ESPRESSO.yml"
    script:


rule ESPRESSO_Q:
    input:
    output:
    log:
        "logs/ESPRESSO_Q.log"
    conda:
        "../envs/ESPRESSO.yml"
    script:
