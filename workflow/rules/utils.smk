import os
import re


rule dump_versions:
    output:
        ver = "versions.txt"
    #TODO: check whether code can be reverted to use workflow.source_path
    #conda: workflow.source_path("envs/env.yml")
    conda: "../envs/env.yml"
    shell:"""
    conda list > {output.ver} 
    """

rule info: ## print pipeline information
    params:
        name = config["workflow"],
        wdir = os.getcwd(),
        repo = config["repo"],
        res  = config["resdir"],
    run:
        print("Pipeline name: ", params.name)
        print("Pipeline working directory: ", params.wdir)
        print("Pipeline results directory: ", params.res)
        print("Pipeline repository: ", params.repo)
    
