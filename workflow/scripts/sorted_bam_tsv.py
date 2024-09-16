import pandas as pd
import sys
import os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sorted_bams = snakemake.input.sorted_bams
samples = snakemake.params.samples

abs_path=[os.path.abspath(bam) for bam in sorted_bams]

df = pd.DataFrame({"sorted_bam_path": abs_path, "sample_name": samples})

df.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
