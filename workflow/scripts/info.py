import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Pipeline name: ", params.name)
print("Pipeline working directory: ", params.wdir)
print("Pipeline repository: ", params.repo)