rule de_analysis:
    input:
        all_counts=rules.merge_counts.output,
    output:
        dispersion_graph=report("de_analysis/dispersion_graph.svg", "../report/dispersion_graph.rst"),
        ma_graph=report("de_analysis/ma_graph.svg","../report/ma_graph.rst"),
        de_heatmap=report("de_analysis/heatmap.svg","../report/heatmap.rst"),
        correlation_matrix=report("de_analysis/correlation_matrix.svg","../report/correlation_matrix.rst"),
        normalized_counts="de_analysis/normalized_counts.csv",
        de_top_heatmap=report("de_analysis/heatmap_top.svg","../report/heatmap_top.rst"),
        lfc_analysis="de_analysis/lfc_analysis.csv",
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
