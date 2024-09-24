localrules:
    reads_manifest,


rule reads_manifest:
    output:
        "flair/reads_manifest.tsv",
    params:
        samples=samples,
        inputdir=config["inputdir"],
    log:
        "logs/flair/samples_tab.log",
    conda:
        "../envs/base.yml"
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
    output:
        "flair/beds/{sample}.bed",
    log:
        "logs/flair/bam2bed_{sample}.log",
    conda:
        "../envs/flair.yml"
    shell:
        "bam2Bed12 --iput_bam {input.sbam}"  # --keep_supplementary to keep supplementary alignments?


rule flair_collapse:
    input:
        beds="flair/beds/{sample}.bed",
        transcriptome="transcriptome/transcriptome.fa",
        annotation="references/standardized_genomic.gtf",
        sample="filter/{sample}_filtered.fq",
    output:
        isob="flair/{sample}/collapse/isoforms.bed",
        isog="flair/{sample}/collapse/isoforms.gtf",
        isof="flair/{sample}/collapse/isoforms.fa",
    params:
        outdir=lambda wildcards: f"flair/{wildcards.sample}/collapse",
    log:
        "logs/flair/{sample}_collapse.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair collapse --genome {input.transcriptome} -gtf {input.annotation}
         --query {input.beds} --reads {input.sample} --output {params.outdir}
         --threads {resources.cpus_per_task}
        """


rule flair_quantify:
    input:
        reads_manifest="flair/reads_manifest.tsv",
        isof=expand("flair/{sample}/collapse/isoforms.fa",sample=samples["sample"]),
    output:
        counts_matrix="flair/quantify/counts_matrix.tsv",
    params:
        outdir=lambda wildcards, output: output[0][:-18],
    log:
        "logs/flair/quantify.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair quantify --reads_manifest {input.reads_manifest}
         --isoforms {input.isof} --output {params.outdir}
         --threads {resources.cpus_per_task}
        """


rule flair_diffexp:
    input:
        counts_matrix="flair/{sample}/quantify/counts_matrix.tsv",
    output:
        genes_deseq2="flair/{sample}/diffexp/genes_deseq2_MCF7_v_A549.tsv",
        genes_deseq2_QCplots="flair/{sample}/diffexp/genes_deseq2_QCplots_MCF7_v_A549.pdf",
        isoforms_deseq2="flair/{sample}/diffexp/isoforms_deseq2_MCF7_v_A549.tsv",
        isoforms_deseq2_QCplots="flair/{sample}/diffexp/isoforms_deseq2_QCplots_MCF7_v_A549.pdf",
        isoforms_drimseq="flair/{sample}/diffexp/isoforms_drimseq_MCF7_v_A549.tsv",
    params:
        outdir=lambda wildcards: f"flair/{wildcards.sample}/diffexp",
    log:
        "logs/flair/{sample}_diffexp.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair_diffexp --counts_matrix {input.counts_matrix}
         --out_dir {params.outdir} --threads {resources.cpus_per_task}
        """


rule flair_diffsplice:
    input:
        isob="flair/{sample}/collapse/isoforms.bed",
        counts_matrix="flair/{sample}/quantify/counts_matrix.tsv",
    output:
        ds_alt3="flair/{sample}/diffsplice/diffsplice.alt3.events.quant.tsv",
        ds_alt5="flair/{sample}/diffsplice/diffsplice.alt5.events.quant.tsv",
        ds_es="flair/{sample}/diffsplice/diffsplice.es.events.quant.tsv",
        ds_ir="flair/{sample}/diffsplice/diffsplice.ir.events.quant.tsv",
        dr_alt3="flair/{sample}/diffsplice/drimseq_alt3_A_v_B.tsv",
        dr_alt5="flair/{sample}/diffsplice/drimseq_alt5_A_v_B.tsv",
        dr_es="flair/{sample}/diffsplice/drimseq_es_A_v_B.tsv",
        dr_ir="flair/{sample}/diffsplice/drimseq_ir_A_v_B.tsv",
    params:
        outdir=lambda wildcards: f"flair/{wildcards.sample}/diffsplice",
    log:
        "logs/flair/{sample}_diffsplice.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair_diffSplice --isoforms {input.isob} --counts_matrix {input.counts_matrix} 
         --out_dir {params.outdir} --threads {resources.cpus_per_task}
        """


# rule plot_isoforms:
#    input:
#        isob="flair/{sample}/collapse/isoforms.bed",
#        counts_matrix="flair/{sample}/quantify/counts_matrix.tsv",
#        gene_name="",
#    output:
#        ""
#    params:
#
#    log:
#        "logs/flair/plot_isoforms.log"
#    conda:
#        "../envs/flair.yml"
#    shell:
#        """
#        plot_isoform_usage {input.isob} {input.counts_matrix} {gene_name}
#         --out_dir {params.outdir} --threads {resources.cpus_per_task}
#        """
