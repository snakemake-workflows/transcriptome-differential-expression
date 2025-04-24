#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pysam
import re
from collections import defaultdict

sys.stderr = sys.stdout = open(snakemake.log[0], "wt")

isoforms_bed = snakemake.input.isob
sorted_sam = snakemake.input.ssam
output_file = snakemake.output.sample_counts
quality_threshold = snakemake.params.qscore

num_match_in_ss_window = snakemake.config["isoform_analysis"]["quantify"][
    "num_match_in_ss_window"
]
large_indel_tolerance = snakemake.config["isoform_analysis"]["quantify"][
    "large_indel_tolerance"
]



def getannotinfo(bedfile):
    transcripttoexons = {}
    with open(bedfile) as bed:
        for line in bed:
            line = line.rstrip().split("\t")
            name, left, right, chrom = line[3], int(line[1]), int(line[2]), line[0]
            blocksizes = [int(n) for n in line[10].rstrip(",").split(",")]
            if line[5] == "+":
                transcripttoexons[name] = blocksizes
            else:
                transcripttoexons[name] = blocksizes[::-1]
    return transcripttoexons


def check_singleexon(read_start, read_end, tlen):
    return read_start < 25 and (read_end > tlen - 25)


def check_exonenddist(blocksize, disttoend, disttoblock):
    if blocksize < 25:
        return disttoend < 5
    else:
        return disttoblock > 25


def check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen):
    left_coverage = check_exonenddist(
        first_blocksize, read_start, first_blocksize - read_start
    )
    right_coverage = check_exonenddist(
        last_blocksize, tlen - read_end, read_end - (tlen - last_blocksize)
    )
    return left_coverage and right_coverage


def check_stringent(coveredpos, exonpos, tlen, blockstarts, blocksizes):
    read_start, read_end = blockstarts[0], blockstarts[-1] + blocksizes[-1]
    first_blocksize, last_blocksize = exonpos[0], exonpos[-1]
    if len(blocksizes) == 1:
        return check_singleexon(read_start, read_end, tlen)
    else:
        return check_firstlastexon(
            first_blocksize, last_blocksize, read_start, read_end, tlen
        )


def check_splicesites(coveredpos, exonpos, tstart, tend):
    currpos = 0
    for i in range(len(exonpos) - 1):
        elen = exonpos[i]
        currpos += elen
        if tstart < currpos < tend:
            ssvals = coveredpos[currpos - 3 : currpos + 3]
            totinsert = sum([x for x in ssvals if x > 1])
            totmatch = sum([x for x in ssvals if x == 1])
            if totmatch - totinsert <= num_match_in_ss_window:
                return False
    return True


def get_matchvals(md):
    matchvals = []
    mdblocks = re.findall(r"\d+|\D+", md)
    for b in mdblocks:
        if b[0] != "^":
            if b.isnumeric():
                matchvals.extend([1] * int(b))
            else:
                matchvals.append(0)
    return matchvals


def process_cigar(matchvals, cigarblocks, startpos):
    matchpos = 0
    coveredpos = [0] * (startpos - 1)
    queryclipping = []
    tendpos = startpos - 1
    blockstarts, blocksizes = [], []
    for btype, blen in cigarblocks:
        if btype in {4, 5}:
            queryclipping.append(blen)
        elif btype == 0:
            coveredpos.extend(matchvals[matchpos : matchpos + blen])
            blockstarts.append(tendpos)
            blocksizes.append(blen)
            matchpos += blen
            tendpos += blen
        elif btype in {2, 3}:
            coveredpos.extend([0] * blen)
            tendpos += blen
            if blen > large_indel_tolerance:
                return True, None, None, None, None, None
        elif btype == 1:
            coveredpos[-1] += blen
            if blen > large_indel_tolerance:
                return True, None, None, None, None, None
    return False, coveredpos, queryclipping, blockstarts, blocksizes, tendpos


def checktranscriptinannot(exondict, tname):
    try:
        exoninfo = exondict[tname]
    except KeyError:
        raise Exception(
            "The transcript names in the annotation fasta do not appear to match the ones in the isoforms file."
        )
    except Exception as ex:
        raise Exception("**check_splice FAILED for %s" % (tname)) from ex
    return exoninfo


def getbesttranscript(tinfo, transcripttoexons):
    passingtranscripts = []
    for tname in tinfo:
        thist = tinfo[tname]
        matchvals = get_matchvals(thist.md)
        indel, coveredpos, queryclip, blockstarts, blocksizes, tendpos = process_cigar(
            matchvals, thist.cigar, thist.startpos
        )
        if not indel:
            exoninfo = checktranscriptinannot(transcripttoexons, tname)
            passesstringent = check_stringent(
                coveredpos, exoninfo, thist.tlen, blockstarts, blocksizes
            )
            if not [exoninfo, passesstringent]:
                continue
            else:
                passingtranscripts.append(
                    [-1 * thist.alignscore, sum(queryclip), thist.tlen, tname]
                )
    if len(passingtranscripts) > 0:
        passingtranscripts.sort()
        return passingtranscripts[0][-1]
    else:
        return None


class IsoAln:
    def __init__(
        self, name=None, q=None, p=None, cigar=None, tlen=None, als=None, md=None
    ):
        self.name = name
        self.q = q
        self.startpos = p
        self.cigar = cigar
        self.tlen = tlen
        self.alignscore = als
        self.md = md


def parsesam(samfile, transcripttoexons):
    lastread = None
    curr_transcripts = {}
    transcripttoreads = defaultdict(list)
    samfile = pysam.AlignmentFile(samfile, "r")
    for read in samfile.fetch(until_eof=True):
        if read.is_mapped:
            readname = read.query_name
            transcript = read.reference_name
            quality = read.mapping_quality
            if quality > quality_threshold:
                pos = read.reference_start
                try:
                    alignscore = read.get_tag("AS")
                    mdtag = read.get_tag("MD")
                except KeyError as ex:
                    raise Exception(
                        f"Missing AS or MD tag in alignment of '{read.query_name}' in '{sorted_sam}'"
                    ) from ex
                cigar = read.cigartuples
                tlen = (
                    samfile.get_reference_length(transcript)
                    if transcript in samfile.references
                    else 0
                )
                if lastread and readname != lastread:
                    assignedt = getbesttranscript(curr_transcripts, transcripttoexons)
                    if assignedt:
                        transcripttoreads[assignedt].append(lastread)
                    curr_transcripts = {}
                curr_transcripts[transcript] = IsoAln(
                    transcript, quality, pos, cigar, tlen, alignscore, mdtag
                )
                lastread = readname
        if lastread:
            assignedt = getbesttranscript(curr_transcripts, transcripttoexons)
            if assignedt:
                transcripttoreads[assignedt].append(lastread)
    samfile.close()
    return transcripttoreads


def write_output(transcripttoreads, outfile):
    with open(outfile, "w") as f:
        for t in transcripttoreads:
            f.write(f"{t}\t{float(len(transcripttoreads[t])):.1f}\n")
    print(f"Output written to: {output_file}")
    for t, reads in transcripttoreads.items():
        print(f"Transcript {t} -> {len(reads)} reads")


transcripttoexons = getannotinfo(isoforms_bed)
transcripttoreads = parsesam(sorted_sam, transcripttoexons)
write_output(transcripttoreads, output_file)
