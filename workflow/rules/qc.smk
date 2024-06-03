import os

configfile: "config/config.yml"
inputdir="/lustre/project/m2_zdvhpc/transcriptome_data/"

# QC and metadata with NanoPlot

if config["summary"] == "None":
    sample_QC=expand("QC/NanoPlot/{sample}.tar.gz", sample=samples["sample"]),
else:
    sample_QC="QC/NanoPlot/summary.tar.gz"


if config["summary"] == "None":

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
            ## max of 39 for our SLURM partition
            cpus_per_task=min(8, 39),  #problem with max(len(input.fastq),39)
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --tsv_stats -f svg "
            "--fastq {input.fastq} -o {output} 2> {log}"

##TODO: string needs to be reworked for variables (aggregate_input in commons.smk)
    rule plot_all_samples:
        input:
            aggregate_input(samples["sample"]),
            #fastq= [os.path.join(config["inputdir"], sample+".fq.gz") for sample in samples["sample"]],
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
            all_samples=rules.plot_all_samples.output,
        output:
            "QC/NanoPlot/{sample}.tar.gz",
        log:
            "logs/NanoPlot/compress_{sample}.log",
        conda:
            None
        shell:
            "tar zcvf {output} {input} 2> {log}"

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
        shell:
            "tar zcvf {output} {input} 2> {log}"
