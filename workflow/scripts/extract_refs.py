import sys
import os
import zipfile
import glob
import shutil

log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file

ref = snakemake.input[0]
output = snakemake.output[0]
accession = snakemake.params.accession

temp_dir = f"temp_{accession}_{'genome' if output.endswith('.fa') else 'annotation'}"

with zipfile.ZipFile(ref, "r") as zf:
    zf.extractall(temp_dir)

ext = "fna" if output.endswith(".fa") else "gff"

for fname in glob.glob(f"{temp_dir}/ncbi_dataset/data/{accession}/*"):
    if fname.endswith(ext):
        shutil.copyfile(fname, output)
        break
else:
    shutil.rmtree(temp_dir)
    raise FileNotFoundError(
        f"No {ext} file found in {ref} for accession no {accession}"
    )

shutil.rmtree(temp_dir)
