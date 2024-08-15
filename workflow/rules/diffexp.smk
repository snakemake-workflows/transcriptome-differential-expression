localrules:
    write_de_params,
    de_analysis,


rule write_de_params:
    output:
        de_params="de_analysis/de_params.tsv",
    log:
        "logs/de_params.log",
    conda:
        "envs/env.yml"
    script:
        "../scripts/de_params.py"


rule de_analysis:
    input:
        de_params="de_analysis/de_params.tsv",
        all_counts=expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"]),
    output:
        dispersion_graph="de_analysis/dispersion_graph.svg",
        ma_graph="de_analysis/ma_graph.svg",
        de_heatmap="de_analysis/heatmap.svg",
        correlation_matrix="de_analysis/correlation_matrix.svg",
        normalized_counts="de_analysis/normalized_counts.csv",
        de_top_heatmap="de_analysis/heatmap_top.svg",
        lfc_analysis="de_analysis/lfc_analysis.csv",
    params:
        samples=samples,
    log:
        "logs/de_analysis.log",
    threads: 4
    conda:
        "envs/env.yml"
    script:
        "../scripts/de_analysis.py"
