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

def generate_help(sfile):
    """Parse out target and help message from file."""
    handler = open(sfile, "r")
    for line in handler:
        match = re.match(r'^rule\s+([a-zA-Z_-]+):.*?## (.*)$$', line)
        if match:
            target, help = match.groups()
            print("%-20s %s" % (target, help))

#rule help: ## print list of all targets with help
#    input:
#        workflow.included
#    run:
#        print("--------------------------------------------")
#        [generate_help(sfile) for sfile in input]
#        print("--------------------------------------------")
rule info: ## print pipeline information
    params:
        name = config["pipeline"],
        wdir = os.getcwd(),
        repo = config["repo"],
        res  = config["resdir"],
    run:
        print("Pipeline name: ", params.name)
        print("Pipeline working directory: ", params.wdir)
        print("Pipeline results directory: ", params.res)
        print("Pipeline repository: ", params.repo)
    
