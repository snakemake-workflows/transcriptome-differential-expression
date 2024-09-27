rule sam_to_bam:
    input:
        sam="alignments/{sample}.sam",
    output:
        "alignments/{sample}.bam",
    log:
        "logs/samtools/samview_{sample}.log",
    params:
        extra=f'{config["sview_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/view"


rule bam_sort:
    input:
        bam="alignments/{sample}.bam",
    output:
        "sorted_alignments/{sample}_sorted.bam",
    log:
        "logs/samtools/samsort_{sample}.log",
    params:
        extra=f'{config["ssort_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/sort"
