import os
from pathlib import Path
import sys
from itertools import product

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError

# global list of valid suffixes
exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

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
    path = Path(os.path.join(config["inputdir"], sample))
    for extension in exts:
        if os.path.exists(path.with_suffix(extension)):
            return path.with_suffix(extension)

    raise WorkflowError(
        f"No valid sample found for sample: '{sample}' with possible extension '{exts}'"
    )


def aggregate_input(samples):
    # possible extensions:
    valids = list()
    for sample, ext in product(samples, exts):
        path = Path(os.path.join(config["inputdir"], sample))

        if os.path.exists(path.with_suffix(ext)):
            valids.append(path.with_suffix(ext))

    if not len(valids):
        raise WorkflowError(f"no valid samples found, allowed extensions are: '{exts}'")
    return valids


def rule_all_input():
    all_input = list()
    all_input.append("versions.txt")
    all_input.extend(expand("QC/NanoPlot/{sample}.tar.gz", sample=samples["sample"]))
    all_input.append("QC/NanoPlot/all_samples.tar.gz")
    all_input.extend(expand("QC/samstats/{sample}.txt", sample=samples["sample"]))
    all_input.extend(expand("QC/qualimap/{sample}.tar.gz", sample=samples["sample"]))
    all_input.extend(
        expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"])
    )
    all_input.append("merged/all_counts.tsv")
    all_input.append("de_analysis/de_params.tsv")
    all_input.append("de_analysis/dispersion_graph.svg")
    all_input.append("de_analysis/ma_graph.svg")
    all_input.append("de_analysis/heatmap.svg")
    all_input.append("de_analysis/lfc_analysis.csv")
    return all_input
