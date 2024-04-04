import glob
import os

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

# samples = (
#    pd.read_csv(workflow.source_path(config["samples"]),
#                                     dtype={"sample": str},
#                                     header=0,
#                                     comment="#")
#    .set_index("sample", drop=False)
#    .sort_index()
# )


def get_fastq_inputs():
    samples, condition, condition2, batch_effect = list(), list(), list(), list()

    inputdir = config["inputdir"]
    extensions = ("fq", "fq.gz", "fastq", "fastq.gz")

    with open(workflow.source_path(config["samples"])) as samplefile:
        samplefile.readline()  # disregard the header
        for line in samplefile:
            if line.startswith("#"):
                continue
            vals = line.split(",")
            samples.append(vals[0].strip())
            condition.append(vals[1].strip())
            condition2.append(vals[2].strip())
            batch_effect.append(vals[3].strip())

    # get sample files
    sample_files = glob_wildcards(os.path.join(config["inputdir"], "{samples}")).samples
    full_paths = list()
    for fname in sample_files:
        path = os.path.join(inputdir, fname)
        for sample in samples:
            if sample in fname:
                full_paths.append(path)
    # attribution to full file path
    all_samples = {
        sample: sample_file for sample, sample_file in zip(samples, full_paths)
    }

    return all_samples, samples, condition, condition2, batch_effect
