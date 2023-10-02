import os
import sys
import pickle as pkl

import matplotlib
matplotlib.use("Agg") # suppress creating of interactive plots
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

print("ready importing")
ncpus = snakemake["config"]["resources.cpus_per_task"]

counts_df = pd.read_csv('snakemake.input.tsv', sep = '\t', header=0)
# we have a header line containing "Reference" as attribute, hence the following line
# otherwise, we would add an index row, with which we cannot work
counts_df.set_index("Reference", inplace = True)
counts_df = counts_df.T
metadata = pd.read_csv(snakemake.input.coldata, sep = '\t', header=0, index_col=0)

# filtering low quality samples, first those with NAs
#samples_to_keep = ~metadata.condition.isna()
#counts_df = counts_df.loc[samples_to_keep]
#metadata = metadata.loc[samples_to_keep]

#TODO: make this configurable
# next we filter out counts, with counts lower than 10
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors=["condition"],
    #refit_cooks=True,
    n_cpus=ncpus,
)
# fit dispersion(s) and LFCs
dds.deseq2()

dds.plot_dispersions(save_path=snakemake.output.dispersion_graph)

# compute p-values for our dispersion
stat_res = DeseqStats(dds, n_cpus=ncpus)

# print the summary
stat_res.summary()

#TODO save to file, this is the template code:
#stat_res.results_df.to_csv(os.path.join(OUTPUT_PATH, "results.csv"))

# performing LFC shrinkage
# a_condition:
a_condition = snakemake.config["condition_a_identifier"]
b_condition = snakemake.config["condition_b_identifier"]
stat_res.lfc_shrink(coeff=f"condition_{b_condition}_vs_{a_condition}")

stat_res.summary()
#TODO: make graph a snakemake target
stat_res.plot_MA(s=20, save_path=snakemake.output.ma_graph)


stat_df = stat_res.results_df

#sns.heatmap(stat_df, cmap='RdYlGn_r', annot=True)
# we need to have a list of all samples:

sns.clustermap(stat_df, xticklables = snakemake.config[""])
plt.savefig(snakemake.output.de_heatmap)
