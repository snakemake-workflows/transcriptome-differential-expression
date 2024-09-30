import pandas as pd
import sys
import os
from pathlib import Path

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

samples_df=snakemake.params.samples[["sample", "condition", "batch"]]
inputdir=snakemake.params.inputdir

exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

def get_sample_path(sample_name, inputdir, exts):
    path=Path(os.path.join(inputdir,sample_name))
    for ext in exts:
        if os.path.exists(str(path)+ext):
            return str(path)+ext
    return "File not found"

samples_df["sample_path"]=samples_df["sample"].apply(lambda x: get_sample_path(x, snakemake.params.inputdir, exts))

samples_df.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
