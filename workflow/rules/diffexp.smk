rule de_analysis:
    input:
        all_counts=rules.merge_counts.output,
    output:
        dispersion_graph=report(
            "de_analysis/dispersion_graph.svg",
            category="Results",
            caption="../report/dispersion_graph.rst",
            labels={
                "figure": "Dispersion graph",
            },
        ),
        ma_graph=report(
            "de_analysis/ma_graph.svg",
            category="Results",
            caption="../report/ma_graph.rst",
            labels={
                "figure": "Bland-Altman plot",
            },
        ),
        de_heatmap=report(
            "de_analysis/heatmap.svg",
            category="Results",
            caption="../report/heatmap.rst",
            labels={
                "figure": "Heatmap",
            },
        ),
        correlation_matrix=report(
            "de_analysis/correlation_matrix.svg",
            category="Results",
            caption="../report/correlation_matrix.rst",
            labels={
                "figure": "Correlation matrix",
            },
        ),
        normalized_counts="de_analysis/normalized_counts.csv",
        de_top_heatmap=report(
            "de_analysis/heatmap_top.svg",
            category="Results",
            caption="../report/heatmap_top.rst",
            labels={
                "figure": "Top heatmap",
            },
        ),
        lfc_analysis="de_analysis/lfc_analysis.csv",
        volcano_plot=report(
            "de_analysis/volcano_plot.svg",
            category="Results",
            caption="../report/volcano_plot.rst",
            labels={
                "figure": "Volcano plot",
            },
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
