 Volcano plot of all genes with their relative gene expression as log2FoldChange (x-axis) and their significance with adjusted p-values (y-axis). 
 The expression strength criteria (dotted lines) lfc_null = {{ snakemake.config["deseq2"]["lfc_null"] }} and the significance threshold (grey line) alpha = {{ snakemake.config["deseq2"]["alpha"] }} determine differentially expressed genes.
 Genes that exceed both are coloured green for overexpressed genes and red for underexpressed genes, other genes are not considered to be differentially expressed between conditions and are coloured grey.
