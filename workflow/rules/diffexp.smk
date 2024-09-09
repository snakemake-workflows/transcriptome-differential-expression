rule de_analysis:
    input:
        all_counts=rules.merge_counts.output,
    output:
        dispersion_graph="de_analysis/dispersion_graph.svg",
        ma_graph="de_analysis/ma_graph.svg",
        de_heatmap="de_analysis/heatmap.svg",
        correlation_matrix="de_analysis/correlation_matrix.svg",
        normalized_counts="de_analysis/normalized_counts.csv",
        de_top_heatmap="de_analysis/heatmap_top.svg",
        lfc_analysis="de_analysis/lfc_analysis.csv",
        de_genes="de_analysis/de_genes.csv",
        volcano_plot="de_analysis/volcano_plot.svg",
    params:
        samples=samples,
    log:
        "logs/de_analysis.log",
    threads: 4
    conda:
        "../envs/pydeseq2.yml"
    script:
        "../scripts/de_analysis.py"
