import os
from glob import glob
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
        dtype={"sample": str, "condition": str, "condition2": str, "batch": str},
        header=0,
        comment="#",
    )
    .set_index("sample", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")


# Verifiy user provided reference files and set input paths
def get_reference_files(config):
    ref = config.get("ref", {})
    accession = ref.get("accession")
    # Validate genome and annotation paths
    genome = ref.get("genome") if ref.get("genome") and os.path.exists(ref["genome"]) and ref["genome"].endswith((".fa", ".fna", ".fasta")) else None
    annotation = ref.get("annotation") if ref.get("annotation") and os.path.exists(ref["annotation"]) and ref["annotation"].endswith((".gtf", ".gff")) else None

    if genome and annotation:
        return {"genome": genome, "annotation": annotation}

    if accession:
        if genome:
            return {"genome": genome}
        if annotation:
            return {"annotation": annotation}
        return {}

    # Throw errors if reference data is invalid
    if genome or annotation:
        raise ValueError("Only one reference file provided, provide either both genome and annotation or an NCBI accession number.")
    raise ValueError("No valid reference files or accession number provided.")


def get_mapped_reads_input(sample):
    path = Path(os.path.join(config["inputdir"], sample))
    for extension in exts:
        # we need to append the extension with +, because
        # path.with_suffix might consider everything after a . in
        # the file name a suffix!
        if os.path.exists(str(path) + extension):
            return str(path) + extension

    raise WorkflowError(
        f"No valid sample found for sample: '{sample}' with possible extension '{exts}'"
    )


def aggregate_input(samples):
    # possible extensions:
    valids = list()
    for sample, ext in product(samples, exts):
        path = Path(os.path.join(config["inputdir"], sample))
        # we need to append the extension with +, because
        # path.with_suffix might consider everything after a . in
        # the file name a suffix!
        if os.path.exists(str(path) + ext):
            valids.append(str(path) + ext)

    if not len(valids):
        raise WorkflowError(f"no valid samples found, allowed extensions are: '{exts}'")
    return valids


def rule_all_input():
    all_input = list()
    all_input.append("versions.txt")
    all_input.extend(
        expand("NanoPlot/{sample}/NanoPlot-report.html", sample=samples["sample"])
    )
    all_input.append("NanoPlot/all_samples/NanoPlot-report.html")
    all_input.extend(expand("QC/bamstats/{sample}.txt", sample=samples["sample"]))
    all_input.extend(
        expand("qualimap/{sample}/qualimapReport.html", sample=samples["sample"])
    )
    all_input.extend(
        expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"])
    )
    all_input.append("merged/all_counts.tsv")
    all_input.append(f"de_analysis/dispersion_graph.{config['deseq2']['figtype']}")
    all_input.append(f"de_analysis/ma_graph.{config['deseq2']['figtype']}")
    all_input.append(f"de_analysis/heatmap.{config['deseq2']['figtype']}")
    all_input.append("de_analysis/lfc_analysis.csv")
    return all_input
