localrules:
    get_references,
    get_annotation,
    get_genome,


rule get_references:
    output:
        # we need two different output, to ensure simultaneous access
        # by the two downstream rules:
        temp("references/ncbi_dataset_a.zip"),
        temp("references/ncbi_dataset_b.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_references.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include gff3,genome &> {log} && mv ncbi_dataset.zip {output[0]};
        cp {output[0]} {output[1]}
        """


rule get_genome:
    input:
        lambda wildcards: get_reference_files(config).get(
            "genome", "references/ncbi_dataset_a.zip"
        ),
    output:
        temp("references/genomic.fa"),
    priority: 10
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_genome.log",
    # conda:
    #    "../envs/reference.yml"
    run:
        import glob, shutil, zipfile

        # we are dealing with a download based on an accession number
        if input[0] == "references/ncbi_dataset_a.zip":
            with zipfile.ZipFile(input[0]) as zf:
                zf.extractall("ncbi_dataset_a")
            for fname in glob.glob(
                f"ncbi_dataset_a/ncbi_dataset/data/{params.accession}/*"
            ):
                if fname.endswith("fna"):
                    shutil.copyfile(fname, f"{output[0]}")
                    shutil.rmtree("ncbi_dataset_a")
                    break
        else:  # a file has been indicated
            shutil.copyfile(input[0], output)


rule get_annotation:
    input:
        lambda wildcards: get_reference_files(config).get(
            "annotation", "references/ncbi_dataset_b.zip"
        ),
    output:
        temp("references/genomic.gff"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_annotation.log",
    # conda:
    #    "../envs/reference.yml"
    run:
        import glob, shutil, zipfile

        # we are dealing with a download based on an accession number
        if input[0] == "references/ncbi_dataset_b.zip":
            with zipfile.ZipFile(input[0]) as zf:
                zf.extractall("ncbi_dataset_b")
            for fname in glob.glob(
                f"ncbi_dataset_b/ncbi_dataset/data/{params.accession}/*"
            ):
                if fname.endswith("gff"):
                    shutil.copyfile(fname, f"{output[0]}")
                    shutil.rmtree("ncbi_dataset_b")
                    break
        else:  # a file has been indicated
            shutil.copyfile(input[0], output[0])
