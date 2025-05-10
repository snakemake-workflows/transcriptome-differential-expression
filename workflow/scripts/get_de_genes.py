#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sorted_lfc_counts = snakemake.input.sorted_lfc_counts
transcriptome = snakemake.input.transcriptome
gene_list = snakemake.output[0]

df = pd.read_csv(sorted_lfc_counts)
df.drop(
    df[df["log2FoldChange"] <= snakemake.config["deseq2"]["lfc_null"]].index,
    inplace=True,
)

gene_names = set(df["Reference"].str.strip())

filtered_records = []
for record in SeqIO.parse(transcriptome, "fasta"):
    record_id = record.id.split()[0]
    if record_id in gene_names:
        filtered_records.append(record)

with open(gene_list, "w") as out:
    SeqIO.write(filtered_records, out, "fasta")
