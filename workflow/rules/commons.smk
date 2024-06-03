import glob
import os

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate
from itertools import product

inputdir="/lustre/project/m2_zdvhpc/transcriptome_data/"

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
    return glob.glob(os.path.join(config["inputdir"], sample) + "*")[0]

def aggregate_input(samples):
    #possible extensions:
    exts = ['fastq', 'fq', 'fastq.gz', 'fq.gz']
    valids = list()
    for sample, ext in product(samples, exts):
        path = os.path.join(inputdir, sample + "." + ext)
        
        if os.path.exists(path):
           valids.append(path)

    #if not len(valids):
    #   raise WorkflowError("errormessage")
    return valids