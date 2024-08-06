import os
import re


rule dump_versions:
    output:
        ver="versions.txt",
    log:
        "logs/utils/dump_ver.log"
    conda:
        "../envs/env.yml"
    # we are using 'ensureconda' because we are unsure which
    # conda flavour is prefered by the user
    shell:
        """
    eval $(ensureconda) list > {output.ver} 2> {log}
    """


rule info:  ## print pipeline information
    params:
        name=config["workflow"],
        wdir=os.getcwd(),
        repo=config["repo"],
    log:
        "logs/utils/info.log"
    run:
        "scripts/info.py"
