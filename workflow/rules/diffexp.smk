rule de_analysis:
    input:
        all_counts=rules.merge_counts.output,
    output:
        dispersion_graph=report(
            "de_analysis/dispersion_graph." + config["deseq2"]["figtype"],
            category="Results",
            caption="../report/dispersion_graph.rst",
        ),
        ma_graph=report(
            "de_analysis/ma_graph." + config["deseq2"]["figtype"],
            category="Results",
            caption="../report/ma_graph.rst",
        ),
        de_heatmap=report(
            "de_analysis/heatmap." + config["deseq2"]["figtype"],
            category="Results",
            caption="../report/heatmap.rst",
        ),
        correlation_matrix=report(
            "de_analysis/correlation_matrix." + config["deseq2"]["figtype"],
            category="Results",
            caption="../report/correlation_matrix.rst",
        ),
        normalized_counts=report("de_analysis/normalized_counts.csv"),
        de_top_heatmap=report(
            "de_analysis/heatmap_top." + config["deseq2"]["figtype"],
            category="Results",
            caption="../report/heatmap_top.rst",
        ),
        sorted_normalized_counts=report("de_analysis/sorted_normalized_counts.csv"),
        lfc_analysis="de_analysis/lfc_analysis.csv",
        volcano_plot=report(
            "de_analysis/volcano_plot." + config["deseq2"]["figtype"],
            category="Results",
            caption="../report/volcano_plot.rst",
        ),
    params:
        samples=samples,
    log:
        "logs/de_analysis.log",
    threads: 4
    conda:
        "../envs/pydeseq2.yml"
    script:
        "../scripts/de_analysis.py"
