rule sam_view:
    input:
        sam="alignments/{sample}.sam",
    output:
        temp("sorted_alignments/{sample}.bam"),
    log:
        "logs/samtools/samview_{sample}.log",
    params:
        extra=f'{config["sview_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/view"


rule sam_sort:
    input:
        sam="alignments/{sample}.sam",
    output:
        temp("sorted_alignments/{sample}_sorted.bam"),
    log:
        "logs/samtools/samsort_{sample}.log",
    params:
        extra=f'{config["ssort_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/sort"
