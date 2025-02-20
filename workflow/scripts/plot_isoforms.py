import sys
import os
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

de_gene_list=snakemake.input.genes[0]
isoforms_bed = snakemake.input.isob
counts_matrix=snakemake.input.counts_matrix
out_dir=snakemake.output[0]

os.makedirs(out_dir, exist_ok=True)

def get_gene_names (de_gene_list):
    df = pd.read_csv(de_gene_list, sep="\t")
    return list(df.iloc[1:,0])

def run_plot_script (isoforms_bed, counts_matrix, gene_name, out_dir):
    cmd = f"plot_isoform_usage {isoforms_bed} {counts_matrix} {gene_name} -o {out_dir}/{gene_name}"
    os.system(cmd)  

for gene in get_gene_names(de_gene_list):
    run_plot_script(isoforms_bed, counts_matrix, gene, out_dir)