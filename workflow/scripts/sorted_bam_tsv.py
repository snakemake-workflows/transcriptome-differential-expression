import pandas as pd
import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sorted_bam_files = snakemake.input.sorted_bams
sample_names = snakemake.params.samples

df = pd.DataFrame({"sorted_bam_path": sorted_bam_files, "sample_name": sample_names})

df.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
