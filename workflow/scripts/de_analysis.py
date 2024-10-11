import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # suppress creating of interactive plots
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from bioinfokit import visuz


from snakemake.exceptions import WorkflowError

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

ncpus = snakemake.threads
samples = snakemake.params.samples

inference = DefaultInference(n_cpus=ncpus)

metadata = samples.loc[:, samples.columns != "samples"]

counts_df = pd.read_csv(f"{snakemake.input.all_counts}", sep="\t", header=0)
# we have a header line containing "Reference" as attribute, hence the following line
# otherwise, we would add an index row, with which we cannot work
counts_df.set_index("Reference", inplace=True)
counts_df = counts_df.T

# next we filter out counts, with counts lower than 10
genes_to_keep = counts_df.columns[
    counts_df.sum(axis=0) >= snakemake.config["deseq2"]["mincount"]
]
counts_df = counts_df[genes_to_keep]

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors=snakemake.config["deseq2"]["design_factors"],
    continuous_factors=snakemake.config["deseq2"].get("continuous_factors", None),
    refit_cooks=True,
    n_cpus=ncpus,
    fit_type=snakemake.config["deseq2"]["fit_type"],
)
# compute normalization factors
# dds.fit_size_factors()
# fit dispersion(s) and LFCs
#  this - fits the size factors
#       - the genewise dispersion
#       - the prior dispersion variance
#       - MAP dispersion
#       - log2fold change
#       - cooks distances
dds.deseq2(fit_type=snakemake.config["deseq2"]["fit_type"])

# fitting genewise dispersions
dds.fit_genewise_dispersions()

dds.plot_dispersions(save_path=f"{snakemake.output.dispersion_graph}")

# compute p-values for our dispersion
stat_res = DeseqStats(dds)
# conditions:
conditions = list()
for condition in samples["condition"]:
    if condition not in conditions:
        conditions.append(condition)

if len(conditions) != 2:
    raise WorkflowError(
        "Only binary conditions are allowed. Make sure your samples.csv only has 2 conditions."
    )
a_condition = conditions[1]  # this order ensures the the lfc shrink condition is met
b_condition = conditions[0]

# run Wald test and plot, perform optional threshold tests, if wanted
if (
    "lfc_null" in snakemake.config["deseq2"]
    or "alt_hypothesis" in snakemake.config["deseq2"]
):
    summary = stat_res.summary(
        lfc_null=snakemake.config["deseq2"].get("lfc_null", 0),
        alt_hypothesis=snakemake.config["deseq2"].get("alt_hypothesis", None),
    )
else:
    summary = stat_res.summary()
# performing LFC shrinkage - we try both combination, because, we
# have no foreknowledge of which conditions comes first
try:
    stat_res.lfc_shrink(coeff=f"condition_{a_condition}_vs_{b_condition}")
except KeyError:
    stat_res.lfc_shrink(coeff=f"condition_{b_condition}_vs_{a_condition}")


stat_res.results_df.to_csv(snakemake.output.lfc_analysis)

stat_res.plot_MA(
    s=snakemake.config["deseq2"]["point_width"],
    save_path=f"{snakemake.output.ma_graph}",
)

# create a clustermap, based on normalized counts
# dds_df = dds.to_df()
# ds_df.to_csv('dds_df.csv'
# getting and applying the scaling factors
sf = dds.obsm["size_factors"]

normalized = counts_df.T * sf

# shorthand for log2fold and pvalue columns
log2foldchange = np.abs(stat_res.results_df["log2FoldChange"])
# 'pvalue' is a pandas series, linear, of length(number of aligned loci)
pvalue = stat_res.results_df["pvalue"]
padj = stat_res.results_df["padj"]

normalized = normalized.join(log2foldchange)
# normalized = normalized.join(np.array(pvalue[1]))
normalized = normalized.join(padj)

# delete rows, which do not meet our p-value criterion
# the comparison operator is >= because we drop all values >= our desired alpha
# normalized.drop(normalized[padj >= snakemake.config["deseq2"]["alpha"]].index, inplace=True)
normalized.to_csv(snakemake.output.normalized_counts)
normalized.to_csv("normalized_my.csv")


# warning: dropping rows before writing might remove the numerical of the selection column
#          the reason is unclear.
sorted = normalized.copy()

sorted.sort_values(by=["log2FoldChange", "padj"], ascending=False, inplace=True)
sorted.drop(sorted[padj >= snakemake.config["deseq2"]["alpha"]].index, inplace=True)
# throw away these columns
# sorted.drop("log2FoldChange", axis=1, inplace=True)
# sorted.drop("padj", axis=1, inplace=True)
sorted.to_csv(snakemake.output.sorted_normalized_counts)

# throw away these columns
normalized.drop("log2FoldChange", axis=1, inplace=True)
normalized.drop("padj", axis=1, inplace=True)
sorted.drop("log2FoldChange", axis=1, inplace=True)
sorted.drop("padj", axis=1, inplace=True)


normalized.dropna(inplace=True)

# precompute linkages, to prevent missing values crashing the script
row_dism = 1 - sorted.T.corr()
row_linkage = hc.linkage(sp.distance.squareform(row_dism), method="complete")
col_dism = 1 - sorted.corr()
col_linkage = hc.linkage(sp.distance.squareform(col_dism), method="complete")

# we calculate the calculation matrix once, with filling NaNs as 0
correlation_matrix = sorted.corr().fillna(0)
# next, a triangle mask is needed, to avoid redundant plotting of the matrix
# in a square
mask = np.triu(np.ones_like(correlation_matrix))

# TODO: add condition labels (e.g. male/female to the map)
cluster = sns.clustermap(
    correlation_matrix,
    cmap=snakemake.config["deseq2"]["colormap"],
    linewidths=0,
)  # , xticklables = metadata.index.to_list())#, yticklabels = sta)
values = cluster.ax_heatmap.collections[0].get_array().reshape(correlation_matrix.shape)
new_values = np.ma.array(values, mask=mask)
cluster.ax_heatmap.collections[0].set_array(new_values)
cluster.ax_col_dendrogram.set_visible(False)
plt.savefig(snakemake.output.correlation_matrix)

# TODO: add condition labels (e.g. male/female to the map)
sns.clustermap(
    normalized.fillna(0),
    cmap=snakemake.config["deseq2"]["colormap"],
    linewidths=0,
    norm=LogNorm(),
)
plt.savefig(snakemake.output.de_heatmap)

n = snakemake.config["deseq2"]["threshold_plot"]
sns.clustermap(
    sorted.fillna(0).iloc[:n],
    cmap=snakemake.config["deseq2"]["colormap"],
    linewidths=0,
    norm=LogNorm(),
)
plt.savefig(snakemake.output.de_top_heatmap)

# our test case has no significant values
# in our CI test, we have no significant data, hence:
if snakemake.config["deseq2"]["alpha"] < 0.9:
    visuz.GeneExpression.volcano(
        df=stat_res.results_df.fillna(1),
        lfc="log2FoldChange",
        pv="padj",
        lfc_thr=(
            snakemake.config["deseq2"]["lfc_null"],
            snakemake.config["deseq2"]["lfc_null"],
        ),
        pv_thr=(
            snakemake.config["deseq2"]["alpha"],
            snakemake.config["deseq2"]["alpha"],
        ),
        sign_line=True,
        gstyle=2,
        show=False,
        plotlegend=True,
        legendpos="upper right",
        legendanchor=(1.46, 1),
        figtype=snakemake.config["deseq2"]["figtype"],
    )
    os.rename(
        "volcano." + snakemake.config["deseq2"]["figtype"],
        snakemake.output.volcano_plot,
    )
else:
    Path(snakemake.output.volcano_plot).touch()
