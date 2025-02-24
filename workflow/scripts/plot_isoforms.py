import sys
import os
import pandas as pd
import subprocess

log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file

de_gene_list = snakemake.input.genes[0]
isoforms_bed = snakemake.input.isob
counts_matrix = snakemake.input.counts_matrix
out_dir = snakemake.output[0]

os.makedirs(out_dir, exist_ok=True)


def get_gene_names(de_gene_list):
    try:
        df = pd.read_csv(de_gene_list, sep="\t")
        if df.empty:
            raise ValueError("Empty gene list file")
        return (gene for gene in df.iloc[:,0])
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        raise ValueError(f"Failed to parse gene list file: {e}")
    except FileNotFoundError:
        raise FileNotFoundError(f"Gene list file not found: {de_gene_list}")


def run_plot_script(isoforms_bed, counts_matrix, gene_name, out_dir):
    try:
        result = subprocess.run(
            [
                "plot_isoform_usage",
                isoforms_bed,
                counts_matrix,
                gene_name,
                "-o",
                f"{out_dir}/{gene_name}",
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Failed to run plot_isoform_usage for gene {gene_name}: {e}"
        )


for gene in get_gene_names(de_gene_list):
    run_plot_script(isoforms_bed, counts_matrix, gene, out_dir)
