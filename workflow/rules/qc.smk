import os


localrules:
    compress_nplot,
    compress_nplot_all,
    compress_map_qc,


configfile: "config/config.yml"


# QC and metadata with NanoPlot
rule plot_samples:
    input:
        fastq=lambda wildcards: get_mapped_reads_input(
            samples["sample"][wildcards.sample]
        ),
    output:
        scatter=report(
            "NanoPlot/{sample}/LengthvsQualityScatterPlot_kde.png",
            category="Quality control",
            caption="../report/nanoplot_samples_scatter_kde.rst",
        ),
    params:
        outdir=directory("NanoPlot/{sample}"),
    log:
        "logs/NanoPlot/{sample}.log",
    resources:
        cpus_per_task=min(len({input}), config["max_cpus"]),  #problem with max(len(input.fastq),39)
    conda:
        "../envs/env.yml"
    shell:
        "NanoPlot --threads {resources.cpus_per_task} --tsv_stats --format png "
        "--fastq {input.fastq} --outdir {params.outdir} 2> {log}"


rule plot_all_samples:
    input:
        aggregate_input(samples["sample"]),
    output:
        scatter=report(
            "NanoPlot/all_samples/LengthvsQualityScatterPlot_kde.png",
            category="Quality control",
            caption="../report/nanoplot_all_scatter_kde.rst",
        ),
    params:
        outdir=directory("NanoPlot/all_samples"),
    log:
        "logs/NanoPlot/all_samples.log",
    conda:
        "../envs/env.yml"
    shell:
        "NanoPlot --threads {resources.cpus_per_task} --tsv_stats --format png "
        "--fastq {input} --outdir {params.outdir} 2> {log}"


rule compress_nplot:
    input:
        samples=rules.plot_samples.output,
    output:
        "QC/NanoPlot/{sample}.tar.gz",
    log:
        "logs/NanoPlot/compress_{sample}.log",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/make_archive.py"


rule compress_nplot_all:
    input:
        all_samples=rules.plot_all_samples.output,
    output:
        "QC/NanoPlot/all_samples.tar.gz",
    log:
        "logs/NanoPlot/compress_all_samples.log",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/make_archive.py"


rule map_qc:
    input:
        bam="sorted_alignments/{sample}_sorted.bam",
    output:
        directory("QC/qualimap/{sample}"),
    log:
        "logs/qualimap/{sample}.log",
    wrapper:
        "v3.13.4/bio/qualimap/bamqc"


rule compress_map_qc:
    input:
        map_qc=rules.map_qc.output,
    output:
        "QC/qualimap/{sample}.tar.gz",
    log:
        "logs/qualimap/compress_{sample}.log",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/make_archive.py"


rule sam_stats:
    input:
        bam="sorted_alignments/{sample}.bam",
    output:
        "QC/samstats/{sample}.txt",
    log:
        "logs/samtools/samstats_{sample}.log",
    params:
        extra=f'{config["sstats_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/stats"
