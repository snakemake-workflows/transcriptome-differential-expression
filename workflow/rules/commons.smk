import fnmatch

def get_fastq_input(inputdir):
    """
    FASTQ-files may end on 'fq', 'fastq' and have an
    additional suffix indicating a compression format
    """
    return [fname for fname in os.listdir(inputdir) if \
            fnmatch.fnmatch(fname, "*[\.fq|\.fastq]*")]