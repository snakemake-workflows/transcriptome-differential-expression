import os


localrules:
    qm_report,


configfile: "config/config.yml"


# QC and metadata with NanoPlot
rule plot_samples:
    input:
        fastq=lambda wildcards: get_mapped_reads_input(
            samples["sample"][wildcards.sample]
        ),
    output:
        scatter=report(
            "NanoPlot/{sample}/NanoPlot-report.html",
            category="Quality control",
            subcategory="NanoPlot",
            caption="../report/nanoplot_sample_report.rst",
            labels={
                "model": "NanoPlot",
                "figure": "{sample}",
            },
        ),
    params:
        outdir=lambda wildcards: f"NanoPlot/{wildcards.sample}",
    log:
        "logs/NanoPlot/{sample}.log",
    resources:
        cpus_per_task=min(len({input}), config["max_cpus"]),  #problem with max(len(input.fastq),39)
    conda:
        "../envs/nanoplot.yml"
    shell:
        "NanoPlot --threads {resources.cpus_per_task} --tsv_stats --format svg "
        "--fastq {input.fastq} --outdir {params.outdir} 2> {log}"


rule plot_all_samples:
    input:
        aggregate_input(samples["sample"]),
    output:
        scatter=report(
            "NanoPlot/all_samples/NanoPlot-report.html",
            category="Quality control",
            subcategory="NanoPlot",
            caption="../report/nanoplot_all_samples_report.rst",
            labels={
                "model": "NanoPlot",
                "figure": "All samples",
            },
        ),
    # This parameter is in line with the Snakemake docs 8.20.3 guideline on how to avoid having parameters as output prefixes
    params:
        outdir=lambda wildcards, output: output[0][:-21],
    log:
        "logs/NanoPlot/all_samples.log",
    conda:
        "../envs/nanoplot.yml"
    shell:
        "NanoPlot --threads {resources.cpus_per_task} --tsv_stats --format svg "
        "--fastq {input} --outdir {params.outdir} 2> {log}"


rule map_qc:
    input:
        bam="sorted_alignments/{sample}_sorted.bam",
    output:
        directory("QC/qualimap/{sample}"),
    log:
        "logs/qualimap/{sample}.log",
    wrapper:
        "v4.4.0/bio/qualimap/bamqc"


# this is a dummy rule to create input for the report because the QualiMap wrapper only accepts directories as valid output
rule qm_report:
    input:
        map_qc=rules.map_qc.output,
    output:
        qm_report=report(
            "qualimap/{sample}/qualimapReport.html",
            category="Quality control",
            caption="../report/qualimap.rst",
        ),
    log:
        "logs/qualimap/{sample}_report.log",
    conda:
        "../envs/base.yml"
    shell:
        "cp -a QC/qualimap/{wildcards.sample} qualimap/ 2> {log}"


rule bam_stats:
    input:
        bam="alignments/{sample}.bam",
    output:
        "QC/bamstats/{sample}.txt",
    log:
        "logs/samtools/bamstats_{sample}.log",
    params:
        extra=config["bamstats_opts"],
    wrapper:
        "v3.13.4/bio/samtools/stats"
