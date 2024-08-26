import sys
from shutil import make_archive
import os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

if isinstance(snakemake.input, str):
    input_dirs = [snakemake.input]
else:
    input_dirs = snakemake.input


base_name = os.path.splitext(os.path.splitext(snakemake.output[0])[0])[0]
root_dir = os.path.dirname(input_dirs[0])
base_dir = os.path.basename(input_dirs[0])
make_archive(base_name, "gztar", root_dir, base_dir)
