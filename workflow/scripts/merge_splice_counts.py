from os import path
import sys
import pandas as pd
import numpy as np
from functools import reduce

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

dfs = {
    x: pd.read_csv(x, sep="\t", header=None, names=["read_id", "Count"])
    for x in snakemake.input
}

ndfs = []
for x, df in dfs.items():
    df["id"] = df["read_id"].apply(lambda s: s.split("_")[-1])
    df["Count"] = df["Count"].astype(float)
    df = df.groupby("id", as_index=False).sum()
    name = path.basename(x).split("_counts.tsv")[0]
    df = df.rename(columns={"Count": name})
    ndfs.append(df[["id", name]])

df_merged = reduce(
    lambda left, right: pd.merge(left, right, on="id", how="outer"), ndfs
)
df_merged = df_merged.fillna(0)
for col in df_merged.columns[1:]:
    df_merged[col] = df_merged[col].astype(float)

df_merged.to_csv(snakemake.output.counts_matrix, sep="\t", index=False)
