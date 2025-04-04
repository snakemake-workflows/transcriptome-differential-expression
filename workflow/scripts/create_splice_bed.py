import os
import sys
import pysam

sys.stderr = sys.stdout = open(snakemake.log[0], "wt")

os.environ["OPENBLAS_NUM_THREADS"] = "1"


def inferMM2JuncStrand(read):
    orientation = read.flag
    try:
        juncDir = read.get_tag("ts")
    except KeyError:
        juncDir = None

    if not juncDir:
        left, right = read.cigar[0], read.cigar[-1]
        s1, s2 = read.seq[:50], read.seq[-50:]
        if ("T" * 10 in s1 and left[0] == 4 and left[1] >= 10) and (
            "A" * 10 in s2 and right[0] == 4 and right[1] >= 10
        ):
            juncDir = "ambig"
        elif "T" * 10 in s1 and left[0] == 4 and left[1] >= 10:
            juncDir = "-" if orientation == 16 else "+"
        elif "A" * 10 in s2 and right[0] == 4 and right[1] >= 10:
            juncDir = "+" if orientation == 16 else "-"
        else:
            juncDir = "ambig"

    else:
        if orientation == 0 and juncDir == "+":
            juncDir = "+"
        elif orientation == 0 and juncDir == "-":
            juncDir = "-"
        elif orientation == 16 and juncDir == "+":
            juncDir = "-"
        elif orientation == 16 and juncDir == "-":
            juncDir = "+"

    return juncDir


def bed_from_cigar(
    alignstart,
    is_reverse,
    cigartuples,
    readname,
    referencename,
    qualscore,
    juncDirection,
):
    positiveTxn = "27,158,119"
    negativeTxn = "217,95,2"
    unknownTxn = "99,99,99"
    refpos = alignstart
    intronblocks = []
    hasmatch = False

    for block in cigartuples:
        if block[0] == 3 and hasmatch:
            intronblocks.append([refpos, refpos + block[1]])
            refpos += block[1]
        elif block[0] in {0, 7, 8, 2}:
            refpos += block[1]
            if block[0] in {0, 7, 8}:
                hasmatch = True

    esizes, estarts = [], [0]
    for i in intronblocks:
        esizes.append(i[0] - (alignstart + estarts[-1]))
        estarts.append(i[1] - alignstart)
    esizes.append(refpos - (alignstart + estarts[-1]))

    rgbcolor = {"+": positiveTxn, "-": negativeTxn, "ambig": unknownTxn}.get(
        juncDirection, unknownTxn
    )

    return [
        referencename,
        str(alignstart),
        str(refpos),
        readname,
        str(qualscore),
        juncDirection,
        str(alignstart),
        str(refpos),
        rgbcolor,
        str(len(intronblocks) + 1),
        ",".join([str(x) for x in esizes]) + ",",
        ",".join([str(x) for x in estarts]) + ",",
    ]


def filter_bam(input_bam, output_bam, output_bed, output_sup, output_trash):
    with pysam.AlignmentFile(input_bam, "rb") as samfile, pysam.AlignmentFile(
        output_bam, "wb", template=samfile
    ) as filtered_bam, pysam.AlignmentFile(
        output_sup, "wb", template=samfile
    ) as sup_bam, pysam.AlignmentFile(
        output_trash, "wb", template=samfile
    ) as trash_bam, open(
        output_bed, "w"
    ) as bed_out:
        for read in samfile.fetch():
            if not read.is_mapped or read.is_secondary:
                trash_bam.write(read)
                continue

            if read.has_tag("SA") or read.is_supplementary:
                sup_bam.write(read)
                continue

            if read.mapping_quality > 1:
                juncstrand = inferMM2JuncStrand(read)
                bedline = bed_from_cigar(
                    read.reference_start,
                    read.is_reverse,
                    read.cigartuples,
                    read.query_name,
                    read.reference_name,
                    read.mapping_quality,
                    juncstrand,
                )
                bed_out.write("\t".join(bedline) + "\n")
                filtered_bam.write(read)
            else:
                trash_bam.write(read)

    pysam.index(output_bam)
    pysam.index(output_sup)
    pysam.index(output_trash)


filter_bam(
    snakemake.input.bam,
    snakemake.output.flair_bam,
    snakemake.output.flair_bed,
    snakemake.output.flair_sup,
    snakemake.output.flair_trash,
)
