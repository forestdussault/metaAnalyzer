import os
from subprocess import Popen

from config import BBMAP_TOOLS_PATH, \
    FASTQC_PATH


def fastq2fasta_bbduk(fastq_filename):
    """
    Converts .fastq to .fasta. Should be done downstream of quality filtering.
    """
    print('\nConverting %s to fasta via --fastq_filter ...' % fastq_filename)

    if '.filtered' in fastq_filename:
        fasta_filename = fastq_filename.replace('.fastq.gz', '.fasta')

        p = Popen('vsearch --fastq_filter {0} -fastaout {1}'.format(fastq_filename, fasta_filename),
                  cwd=os.path.dirname(fastq_filename),
                  shell=True,
                  executable='/bin/bash')
        p.wait()

    else:
        print('No filtered fastq file provided. Skipping %s' % fastq_filename)
        pass


def quality_trim_bbduk(fastq_filename):
    """
    Quality trimming with BBDuk
    Quality score threshold of 10 (Brian Bushnell suggestion for Illumina reads)
    Minimum length of 75 to replicate Pal et al. (2016) though all reads should be ~100bp from MiSeq runs
    """
    fastq_filename_filtered = fastq_filename.replace('.fastq.gz', '.filtered.fastq.gz')
    filepath_in = fastq_filename
    filepath_out = fastq_filename_filtered

    # Quality-trim to Q10 using Phred algorithm with BBDuk. Trims left and right sides of reads.
    p = Popen('./bbduk.sh -Xmx1g in={0} out={1} qtrim=rl trimq=10 '
              'minlen=75 overwrite=true'.format(filepath_in, filepath_out),
              cwd=BBMAP_TOOLS_PATH,
              shell=True,
              executable='/bin/bash')
    p.wait()


def run_fastqc(fastq_filename):
    p = Popen('./fastqc {0} --extract'.format(fastq_filename),
              cwd=FASTQC_PATH,
              shell=True,
              executable='/bin/bash')
    p.wait()