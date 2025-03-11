 The `MA plot, <https://en.wikipedia.org/wiki/MA_plot>`_ is created by computing p-values using Wald-tests.
 The plot compares for each transcript the mean abundance across samples (x-axis) and the log2 foldchange as a ratio of expression between the two conditions (y-axis).
 Genes with significant changes in expression that fall outside of the significance criteria lfc_null = {{ snakemake.config["deseq2"]["lfc_null"] }} are highlighted in red.
