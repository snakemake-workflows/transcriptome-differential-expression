#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
import gzip

# start logging
sys.stderr = sys.stdout = open(snakemake.log[0], "wt")

if snakemake.config["min_length"] == 0:
    snakemake.output = snakemake.input
else:
    with gzip.open(snakemake.input[0], "rt") as sample : 
        input_iterator = SeqIO.parse(sample, "fastq")

        filter_iterator = (read for read in input_iterator if len(read.seq) >= snakemake.config["min_length"])


        SeqIO.write(filter_iterator, snakemake.output[0], "fastq")