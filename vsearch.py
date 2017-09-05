import os
from subprocess import Popen

from config import AMR_DB_PATH, CARB_DB_PATH, DATA_ANALYSIS_PATH

def usearch_global(fasta_filename):
    """
    Run vsearch on a fastq file. No limit to # of hits. Minimum 80% identity.
    """
    query_target = fasta_filename
    ref_db = AMR_DB_PATH
    p = Popen('vsearch --usearch_global {0} --db {1} '
              ' --id 0.8 --maxaccepts 0 --samheader --samout {2}'.format(query_target, CARB_DB_PATH, fasta_filename.replace('.fasta','.sam')),
              shell=True,
              executable="/bin/bash")
    p.wait()


def fastq2fasta(fastq_filename):
    print('Attempting to convert {} to FASTA format...'.format(fastq_filename))
    p = Popen('fastq_to_fasta -v -i {0} -o {1}'.format(fastq_filename, fastq_filename.replace('.fastq','.fasta')),
              cwd='/home/dussaultf/Applications/fastx-toolkit/bin',
              shell=True,
              executable="/bin/bash")
    p.wait()
    print('Completed conversion.')


id_list = os.listdir(DATA_ANALYSIS_PATH)

fastq2fasta('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/4477873.3/4477873.3.050.upload.filtered.fastq')
usearch_global('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/4477873.3/4477873.3.050.upload.filtered.fasta')