import os
import re


rule dump_versions:
    output:
        ver="versions.txt",
    conda:
        "../envs/env.yml"
    # we are using 'ensureconda' because we are unsure which
    # conda flavour is prefered by the user
    shell:
        """
    eval $(ensureconda) list > {output.ver} 
    """


rule info:  ## print pipeline information
    params:
        name=config["workflow"],
        wdir=os.getcwd(),
        repo=config["repo"],
    run:
        print("Pipeline name: ", params.name)
        print("Pipeline working directory: ", params.wdir)
        print("Pipeline repository: ", params.repo)
