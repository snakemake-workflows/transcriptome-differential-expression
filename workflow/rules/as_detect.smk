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
        sbami="sorted_alignments/{sample}_sorted.bam.bai",
    output:
        "flair/beds/{sample}.bed",
    log:
        "logs/flair/bam2bed_{sample}.log",
    conda:
        "../envs/flair.yml"
    shell:
        "bam2Bed12 --input_bam {input.sbam} > {output} &> {log}"  # --keep_supplementary to keep supplementary alignments?

rule concatenate_beds:
    input:
        expand("flair/beds/{sample}.bed",sample=samples["sample"])
    output:
        "flair/beds/all_samples.bed",
    log:
        "logs/flair/concatenate_beds.log"
    conda:
        "../envs/base.yml"
    shell:
        "cat {input} >> {output}"


rule flair_collapse:
    input:
        beds="flair/beds/all_samples.bed",
        transcriptome="transcriptome/transcriptome.fa",
        annotation="references/standardized_genomic.gtf",
        sample=expand("filter/{sample}_filtered.fq",sample=samples["sample"]),
    output:
        isob="flair/isoforms.bed",
        isog="flair/isoforms.gtf",
        isof="flair/isoforms.fa",
    params:
        outdir="flair",
    log:
        "logs/flair/collapse.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair collapse --genome {input.transcriptome} -gtf {input.annotation}
         --query {input.beds} --reads {input.sample} --output {params.outdir}
         --threads {resources.cpus_per_task} &> {log}
        """


rule flair_quantify:
    input:
        reads_manifest="flair/reads_manifest.tsv",
        isof="flair/isoforms.fa",
    output:
        counts_matrix="flair/counts_matrix.tsv",
    params:
        outdir="flair",
    log:
        "logs/flair/quantify.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair quantify --reads_manifest {input.reads_manifest}
         --isoforms {input.isof} --output {params.outdir}
         --threads {resources.cpus_per_task} &> {log}
        """


rule flair_diffexp:
    input:
        counts_matrix="flair/counts_matrix.tsv",
    output:
        genes_deseq2="flair/genes_deseq2_MCF7_v_A549.tsv",
        genes_deseq2_QCplots="flair/genes_deseq2_QCplots_MCF7_v_A549.pdf",
        isoforms_deseq2="flair/isoforms_deseq2_MCF7_v_A549.tsv",
        isoforms_deseq2_QCplots="flair/isoforms_deseq2_QCplots_MCF7_v_A549.pdf",
        isoforms_drimseq="flair/isoforms_drimseq_MCF7_v_A549.tsv",
    params:
        outdir="flair",
    log:
        "logs/flair/diffexp.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair_diffexp --counts_matrix {input.counts_matrix}
         --out_dir {params.outdir} --threads {resources.cpus_per_task}
         &> {log}
        """


rule flair_diffsplice:
    input:
        isob="flair/isoforms.bed",
        counts_matrix="flair/counts_matrix.tsv",
    output:
        ds_alt3="flair/diffsplice.alt3.events.quant.tsv",
        ds_alt5="flair/diffsplice.alt5.events.quant.tsv",
        ds_es="flair/diffsplice.es.events.quant.tsv",
        ds_ir="flair/diffsplice.ir.events.quant.tsv",
        dr_alt3="flair/drimseq_alt3_A_v_B.tsv",
        dr_alt5="flair/drimseq_alt5_A_v_B.tsv",
        dr_es="flair/drimseq_es_A_v_B.tsv",
        dr_ir="flair/drimseq_ir_A_v_B.tsv",
    params:
        outdir="flair",
    log:
        "logs/flair/diffsplice.log",
    conda:
        "../envs/flair.yml"
    shell:
        """
        flair_diffSplice --isoforms {input.isob} --counts_matrix {input.counts_matrix} 
         --out_dir {params.outdir} --threads {resources.cpus_per_task} &> {log}
        """


# rule plot_isoforms:
#    input:
#        isob="flair/isoforms.bed",
#        counts_matrix="flair/counts_matrix.tsv",
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
