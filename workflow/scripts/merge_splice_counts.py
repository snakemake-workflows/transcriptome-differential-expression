from os import path
import sys
import pandas as pd
import numpy as np
from functools import reduce

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

dfs = {x: pd.read_csv(x, sep="\t", header=None, names=["ids", "Count"]) for x in snakemake.input.sample_counts}

ndfs = []
for x, df in dfs.items():
    df.Count = np.array(df.Count, dtype=int)
    df = df[df.Count > 0]
    df = df.sort_values(by=["Count"], ascending=False)

    name = path.dirname(x).rsplit("/", 1)[-1]
    if name.endswith("_counts"):
        name = name.replace("_counts", "")

    df = df.rename(columns={"Count": name})
    df = df[["ids", name]]
    ndfs.append(df)

df_merged = reduce(
    lambda left, right: pd.merge(left, right, on="Reference", how="outer"), ndfs
)

df_merged = df_merged.fillna(0).astype({col: int for col in df_merged.columns if col != "Reference"})
df_merged.to_csv(snakemake.output.samples_counts, sep="\t", index=False)
