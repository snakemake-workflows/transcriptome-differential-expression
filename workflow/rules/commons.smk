import os
from pathlib import Path
import sys
from itertools import product

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(
        os.path.join(
            os.path.realpath(os.path.dirname(workflow.configfiles[0])),
            config["samples"],
        ),
        sep=r"\s+",
        dtype={"sample": str, "condition": str, "condition2": str, "batch_effect": str},
        header=0,
        comment="#",
    )
    .set_index("sample", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")


def get_mapped_reads_input(sample):
    return list(Path(os.path.join(config["inputdir"], sample).glob("*")))[0]


def aggregate_input(samples):
    # possible extensions:
    exts = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    valids = list()
    for sample, ext in product(samples, exts):
        path = Path(os.path.join(config["inputdir"], sample))

        if os.path.exists(path.with_suffix(ext)):
            valids.append(path)

    if not len(valids):
        raise WorkflowError(f"no valid samples found, allowed extensions are: {exts}")
    return valids
