localrules:
    write_coldata,
    write_de_params,
    de_analysis,


rule write_coldata:
    output:
        coldata="de_analysis/coldata.tsv",
    run:
        with open(f"{output}", "w") as outfile:
            outstring = "\t".join(samples.head())
            outfile.write(outstring)


rule write_de_params:
    output:
        de_params="de_analysis/de_params.tsv",
    run:
        d = OrderedDict()
        d["Annotation"] = [config["annotation"]]
        d["min_samps_gene_expr"] = [config["min_samps_gene_expr"]]
        d["min_samps_feature_expr"] = [config["min_samps_feature_expr"]]
        d["min_gene_expr"] = [config["min_gene_expr"]]
        d["min_feature_expr"] = [config["min_feature_expr"]]
        df = pd.DataFrame(d)
        df.to_csv(output.de_params, sep="\t", index=False)


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
