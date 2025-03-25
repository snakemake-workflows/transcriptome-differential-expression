import os
import re


localrules:
    dump_versions,


rule dump_versions:
    output:
        ver="versions.txt",
    log:
        "logs/utils/dump_ver.log",
    conda:
        "../envs/ensureconda.yml"
    # we are using 'ensureconda' because we are unsure which
    # conda flavour is prefered by the user
    shell:
        """
    eval $(ensureconda) list > {output.ver} 2> {log}
    """
