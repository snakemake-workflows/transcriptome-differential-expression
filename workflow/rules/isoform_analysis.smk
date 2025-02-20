localrules:
    reads_manifest,
    concatenate_beds,


# Construct a flair readable TSV file for samples
rule reads_manifest:
    output:
        "iso_analysis/reads_manifest.tsv",
    params:
        samples=samples,
        inputdir=config["inputdir"],
    log:
        "logs/flair/reads_manifest.log",
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/reads_manifest.py"


rule gff_to_gtf:
    input:
        "references/standardized_genomic.gff",
    output:
        "references/standardized_genomic.gtf",
    log:
        "logs/gffread/gff_to_gtf.log",
    conda:
        "../envs/gffread.yml"
    shell:
        """
        gffread -t {resources.cpus_per_task} -E {input} -T -o {output} &> {log}    
        """


rule bam_to_bed:
    input:
        sbam="sorted_alignments/{sample}_sorted.bam",
        sbami="sorted_alignments/{sample}_sorted.bam.bai",
    output:
        "iso_analysis/beds/{sample}.bed",
    log:
        "logs/flair/bam2bed_{sample}.log",
    conda:
        "../envs/flair.yml"
    shell:
        "bam2Bed12 --input_bam {input.sbam} > {output} 2> {log}"  # --keep_supplementary to keep supplementary alignments?


rule concatenate_beds:
    input:
        expand("iso_analysis/beds/{sample}.bed", sample=samples["sample"]),
    output:
        "iso_analysis/beds/all_samples.bed",
    log:
        "logs/flair/concatenate_beds.log",
    conda:
        "../envs/base.yml"
    shell:
        "rm -f {output}; cat {input} >> {output}"


rule build_flair_genome_index:
    input:
        target="references/genomic.fa",
    output:
        index="index/flair_genome_index.mmi",
    params:
        extra=config["minimap2"]["index_opts"],
    log:
        "logs/flair/index.log",
    threads: 4
    wrapper:
        "v3.13.4/bio/minimap2/index"


rule flair_align:
    input:
        genome="references/genomic.fa",
        sample=expand("filter/{sample}_filtered.fq", sample=samples["sample"]),
        index="index/flair_genome_index.mmi",
    output:
        flair_beds="iso_analysis/align/flair.bed",
    params:
        outdir=lambda wildcards, output: output[0][:-4],
    log:
        "logs/flair/align.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair align --reads {input.sample} --genome {input.genome}  \
        --mm_index {input.index} --output {params.outdir} \
        --threads {resources.cpus_per_task} &> {log}
        """


rule flair_correct:
    input:
        genome="references/genomic.fa",
        flair_beds="iso_analysis/align/flair.bed",
        annotation="references/standardized_genomic.gtf",
    output:
        beds_cor="iso_analysis/align/flair_all_corrected.bed",
    params:
        outdir=lambda wildcards, output: output[0][:-18],
    log:
        "logs/flair/correct.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair correct --query {input.flair_beds} --genome {input.genome}  \
        --gtf {input.annotation} --output {params.outdir} \
        --threads {resources.cpus_per_task} &> {log}
        """


rule flair_collapse:
    input:
        beds="iso_analysis/align/flair_all_corrected.bed",
        genome="references/genomic.fa",
        annotation="references/standardized_genomic.gtf",
        sample=expand("filter/{sample}_filtered.fq", sample=samples["sample"]),
    output:
        isob="iso_analysis/collapse/flair.isoforms.bed",
        isof="iso_analysis/collapse/flair.isoforms.fa",
    params:
        outdir=lambda wildcards, output: output[0][:-13],
        qscore=config["FLAIR"]["qscore"],
        opts=config["FLAIR"]["col_opts"],
    log:
        "logs/flair/collapse.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair collapse --genome {input.genome} --gtf {input.annotation} --query {input.beds} \
        --reads {input.sample} --output {params.outdir} --quality {params.qscore} --no_gtf_end_adjustment \
        {params.opts} --threads {resources.cpus_per_task} &> {log}
        """


rule flair_quantify:
    input:
        reads_manifest="iso_analysis/reads_manifest.tsv",
        isof="iso_analysis/collapse/flair.isoforms.fa",
        isob="iso_analysis/collapse/flair.isoforms.bed",
    output:
        counts_matrix="iso_analysis/quantify/flair.counts.tsv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        tmp_dir="iso_analysis/quantify/tmp",
        qscore=config["FLAIR"]["qscore"],
    log:
        "logs/flair/quantify.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair quantify --reads_manifest {input.reads_manifest} --isoforms {input.isof} \
        --isoform_bed {input.isob} --output {params.outdir} --quality {params.qscore} \
        --temp_dir {params.tmp_dir} --stringent --threads {resources.cpus_per_task} \
        &> {log}
        """


## output is not "MCF7" and "A549" but reads_manifest condition "A" vs "B"
rule flair_diffexp:
    input:
        counts_matrix="iso_analysis/quantify/flair.counts.tsv",
    output:
        genes_deseq2="iso_analysis/diffexp/genes_deseq2_{condition_value1}_v_{condition_value2}.tsv",
        genes_deseq2_QCplots="iso_analysis/diffexp/genes_deseq2_QCplots_{condition_value1}_v_{condition_value2}.pdf",
        isoforms_deseq2="iso_analysis/diffexp/isoforms_deseq2_{condition_value1}_v_{condition_value2}.tsv",
        isoforms_deseq2_QCplots="iso_analysis/diffexp/isoforms_deseq2_QCplots_{condition_value1}_v_{condition_value2}.pdf",
        isoforms_drimseq="iso_analysis/diffexp/isoforms_drimseq_{condition_value1}_v_{condition_value2}.tsv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        exp_thresh=config["FLAIR"]["exp_thresh"],
    log:
        "logs/flair/diffexp_{condition_value1}_v_{condition_value2}.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair diffexp --counts_matrix {input.counts_matrix}  --out_dir {params.outdir} \
        --out_dir_force --exp_thresh {params.exp_thresh} --threads {resources.cpus_per_task} \
        &> {log}
        """


rule plot_isoforms:
    input:
        genes=expand(
            "iso_analysis/diffexp/genes_deseq2_{condition_value1}_v_{condition_value2}.tsv",
            condition_value1=condition_value1,
            condition_value2=condition_value2,
        ),
        isob="iso_analysis/collapse/flair.isoforms.bed",
        counts_matrix="iso_analysis/quantify/flair.counts.tsv",
    output:
        out_dir=directory("iso_analysis/plots"),
    log:
        "logs/flair/plot_isoforms.log",
    conda:
        "../envs/flair.yml"
    script:
        "../scripts/plot_isoforms.py"