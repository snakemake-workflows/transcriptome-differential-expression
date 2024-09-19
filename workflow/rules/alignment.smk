rule build_minimap_index:  ## build minimap2 index
    input:
        target="transcriptome/transcriptome.fa",
    output:
        index=temp("index/transcriptome_index.mmi"),
    params:
        extra=config["minimap_index_opts"],
    log:
        "logs/minimap2/index.log",
    wrapper:
        "v3.13.4/bio/minimap2/index"


# mapping reads with minimap2
rule map_reads:
    input:
        target="index/transcriptome_index.mmi",
        query="filter/{sample}_filtered.fq",
    output:
        temp("alignments/{sample}.sam"),
    log:
        "logs/minimap2/mapping_{sample}.log",
    params:
        extra=f"-p {config['secondary_score_ratio']} -N {config['maximum_secondary']} {config['minimap2_opts']}",
    wrapper:
        "v3.13.4/bio/minimap2/aligner"
