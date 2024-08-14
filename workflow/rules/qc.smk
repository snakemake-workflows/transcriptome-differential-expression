import os


localrules:
    compress_nplot,
    compress_nplot_all,
    compress_map_qc,


configfile: "config/config.yml"


# inputdir = config["inputdir"]  # "/lustre/project/m2_zdvhpc/transcriptome_data/"


# QC and metadata with NanoPlot

if not config["summary"].endswith(".txt"):
    sample_QC = (
        (expand("QC/NanoPlot/{sample}.tar.gz", sample=samples["sample"]),),
        "QC/NanoPlot/all_samples.tar.gz",
    )
else:
    sample_QC = "QC/NanoPlot/summary.tar.gz"


if not config["summary"].endswith(".txt"):

    rule plot_samples:
        input:
            fastq=lambda wildcards: get_mapped_reads_input(
                samples["sample"][wildcards.sample]
            ),
        output:
            directory("NanoPlot/{sample}"),
        log:
            "logs/NanoPlot/{sample}.log",
        resources:
            cpus_per_task=min(len({input}), config["max_cpus"]),  #problem with max(len(input.fastq),39)
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --tsv_stats -f svg "
            "--fastq {input.fastq} -o {output} 2> {log}"

    rule plot_all_samples:
        input:
            aggregate_input(samples["sample"]),
        output:
            directory("NanoPlot/all_samples"),
        log:
            "logs/NanoPlot/all_samples.log",
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --tsv_stats -f svg "
            "--fastq {input} -o {output} 2> {log}"

    rule compress_nplot:
        input:
            samples=rules.plot_samples.output,
        output:
            "QC/NanoPlot/{sample}.tar.gz",
        log:
            "logs/NanoPlot/compress_{sample}.log",
        conda:
            "../envs/env.yml"
        shell:
            "tar zcvf {output} {input} &> {log}"

    rule compress_nplot_all:
        input:
            all_samples=rules.plot_all_samples.output,
        output:
            "QC/NanoPlot/all_samples.tar.gz",
        log:
            "logs/NanoPlot/compress_all_samples.log",
        conda:
            "../envs/env.yml"
        shell:
            "tar zcvf {output} {input} &> {log}"

else:

    rule plot_samples:
        input:
            summary=config["summary"],
        output:
            directory("NanoPlot"),
        log:
            "logs/NanoPlot/NanoPlot.log",
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --barcoded --tsv_stats "
            "--summary {input.summary} -o {output} 2> {log}"

    rule compress_nplot:
        input:
            samples=rules.plot_samples.output,
        output:
            "QC/NanoPlot/summary.tar.gz",
        log:
            "logs/NanoPlot/compress_summary.log",
        conda:
            "../envs/env.yml"
        shell:
            "tar zcvf {output} {input} &> {log}"


rule map_qc:
    input:
        bam="sorted_alignments/{sample}_sorted.bam",
    output:
        directory("QC/qualimap/{sample}"),
    log:
        "logs/qualimap/{sample}.log",
    conda:
        "../envs/env.yml"
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
        "../envs/env.yml"
    shell:
        "tar zcvf {output} {input} &> {log}"


rule sam_stats:
    input:
        bam="sorted_alignments/{sample}.bam",
    output:
        "QC/samstats/{sample}.txt",
    log:
        "logs/samtools/samstats_{sample}.log",
    params:
        extra=f'{config["sstats_opts"]}',
    conda:
        "../envs/env.yml"
    wrapper:
        "v3.13.4/bio/samtools/stats"
