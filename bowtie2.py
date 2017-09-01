import os
from subprocess import Popen, check_call

from config import DATABASE_PATH, CARB_DB_PATH, DATA_ANALYSIS_PATH

def make_bowtie2_db(fasta_file, db_name):
    p = Popen('bowtie2-build {0} {1}'.format(fasta_file, db_name),
              cwd=os.path.dirname(fasta_file),
              shell=True,
              executable='/bin/bash')
    p.wait()

# make_bowtie2_db('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/bowtie2/carbapenemase_db/carbapenemase_merged.fasta','carbapenemase_db')

def bowtie2_alignment(database_name, fastq_file, database_path):
    p = Popen('bowtie2 '
              '-x {0} '
              '--very-sensitive-local '
              '-U {1} '
              '-S {2}'.format(database_name, fastq_file, fastq_file.replace('.fastq','.{}.sam'.format(database_name))),
              cwd=database_path,
              shell=True,
              executable='/bin/bash'
              )
    p.wait()

print('\nAligning against amr_db...')
bowtie2_alignment('amr_db',
          '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/4481526.3/4481526.3.050.upload.filtered.fastq',
          '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/bowtie2/amr_db/')

print('\nAligning against carbapenemase_db...')
bowtie2_alignment('carbapenemase_db',
          '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/4481526.3/4481526.3.050.upload.filtered.fastq',
          '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/bowtie2/carbapenemase_db/')