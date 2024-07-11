import os
import sys
import pickle as pkl

import matplotlib

matplotlib.use("Agg")  # suppress creating of interactive plots
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import seaborn as sns
import time

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data



sys.stderr = sys.stdout = open(snakemake.log[0], "w")

ncpus = snakemake.threads
print(snakemake.params.samples)

counts_df = pd.read_csv(f"{snakemake.input.all_counts}", sep="\t", header=0)
# we have a header line containing "Reference" as attribute, hence the following line
# otherwise, we would add an index row, with which we cannot work
counts_df.set_index("sample", inplace=True)
counts_df = counts_df.T
metadata = pd.read_csv(f"{snakemake.input.coldata}", sep=",", header=0, index_col=0)


# TODO: make this configurable
# next we filter out counts, with counts lower than 10
print(counts_df.sum(axis=0))
print(snakemake.config["mincount"])
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= snakemake.config["mincount"]]
print(genes_to_keep)
counts_df = counts_df[genes_to_keep]

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors=["condition"],
    refit_cooks=True,
    n_cpus=ncpus,
)
# fit dispersion(s) and LFCs
#  this - fits the size factors
#       - the genewise dispersion
#       - the prior dispersion variance
#       - MAP dispersion
#       - log2fold change
#       - cooks distances
dds.deseq2()

# compute normalization factors
# dds.fit_size_factors()
# fitting genewise dispersions
# dds.fit_genewise_dispersions()

dds.plot_dispersions(save_path=f"{snakemake.output.dispersion_graph}")

# compute p-values for our dispersion
stat_res = DeseqStats(dds, n_cpus=ncpus)
# a_condition:
a_condition = snakemake.config["condition_a_identifier"]
# b_condition
b_condition = snakemake.config["condition_b_identifier"]

# run Wald test and plot, perform optional threshold tests, if wanted
if "lfc_null" in snakemake.config or "alt_hypothesis" in snakemake.config:
    summary = stat_res.summary(
        lfc_null=snakemake.config.get("lfc_null", 0),
        alt_hypothesis=snakemake.config.get("alt_hypothesis", None),
    )
else:
    summary = stat_res.summary()
# performing LFC shrinkage
stat_res.lfc_shrink(coeff=f"condition_{b_condition}_vs_{a_condition}")

stat_res.results_df.to_csv(snakemake.output.lfc_analysis)

stat_res.plot_MA(
    s=snakemake.config["point_width"], save_path=f"{snakemake.output.ma_graph}"
)

# create a clustermap, based on normalized counts
# dds_df = dds.to_df()
# ds_df.to_csv('dds_df.csv'
# getting and applying the scaling factors
sf = dds.obsm["size_factors"]

normalized = counts_df.values.T * sf

# transpose back
normalized = normalized.T

# append log2fold and pvalue columns
normalized.join(stats_res.results_df["log2FoldChange"])
normalized.join(stats_res.results_df["pvalue"])

# get indices where the matrix is 0 - may happen
# zero_indices = np.argwhere(normalized == 0)
# TODO: try implementing imputation
# for now, we remove those, which are zeroed
normalized = normalized[~np.all(normalized == 0, axis=1)]

# get column names and row names
column_names = counts_df.columns.values.tolist()
row_names = counts_df.index.values.tolist()

normalized = pd.DataFrame(normalized, index=row_names, columns=column_names)

# Order columns according to traits - generally column
# order can be arbitrary, but for the headmap, we want.
a_samples, b_samples = list(), list()
for sample_name in row_names:
    if a_condition in sample_name:
        a_samples.append(sample_name)
    else:
        b_samples.append(sample_name)
assert a_samples, f"list 'a_samples' is empty, '{a_condition}' not unique?"
assert b_samples, f"list 'b_samples' is empty, '{a_condition}' not unique?"
# total list
samples = a_samples + b_samples

# final orientation and order
normalized = normalized.T[samples]


# get the means of our conditions
# a_mean = normalized[a_samples].mean(axis=1)
# b_mean = normalized[b_samples].mean(axis=1)

# a_over_b = a_mean/b_mean
# b_over_a = b_mean/a_mean
# which is bigger?
# ratio = list(None for _ in a_over_b)
# for index, state in enumerate(a_over_b.ge(b_over_a)):
# enter ratio, but check for non-inf-ness, first:
#    print(a_over_b[index], b_over_a[index])
#    if state:
#        if np.isinf(a_over_b[index]):
#            normalized.drop(index, inplace=True)
#        ratio[index] = a_over_b[index]
#    else:
#        if np.isinf(b_over_a[index]):
#            normalized.drop(index, inplace=True)
#       ratio[index] = b_over_a[index]

# normalized["ratio"] = ratio
# normalized = normalized[~np.all(normalized == np.inf, axis=1)]
# print(normalized)
# now sort according to the log2foldchange
normalized.sort_values(by="log2FoldChange")
# delete rows, which do not meet our p-value criterion
normalized.drop(normalized["pvalue"] > 0.05, inplace=True)
# through away this column
# normalized.drop(["ratio"], axis=1, inplace=True)

# normalized.loc[normalized.index.difference(normalized.dropna(how='all').index)]
# print(normalized)

sns.clustermap(
    normalized[samples], cmap=snakemake.config["colormap"], linewidths=0, norm=LogNorm()
)  # , xticklables = metadata.index.to_list())#, yticklabels = sta)
plt.savefig(snakemake.output.de_heatmap)
n = snakemake.config["threshold_plot"]
sns.clustermap(
    normalized.iloc[:n][samples],
    cmap=snakemake.config["colormap"],
    linewidths=0,
    norm=LogNorm(),
)
plt.savefig(snakemake.output.de_top_heatmap)
