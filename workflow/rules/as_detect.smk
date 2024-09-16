localrules:
    samples_tab

rule samples_tab:
    input:
        sorted_bams=expand("sorted_alignments/{sample}_sorted.bam",sample=samples["sample"]),
    output:
        "ESPRESSO/sorted_bam_samples.tsv",
    params:
        samples=samples["sample"],
    log:
        "logs/ESPRESSO/samples_tab.log"
    conda:
        "../envs/base.yml",
    script:
        "../scripts/sorted_bam_tsv.py"
