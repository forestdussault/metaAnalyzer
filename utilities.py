import os
from subprocess import Popen

from config import ROOT_PATH,\
    DATABASE_PATH,\
    ALL_DATASETS,\
    MGRAST_DOWNLOAD_PATH,\
    BBDUK_PATH


def retrieve_id_list(list_of_ids):
    """
    Parses a text file with individual IDs on newlines and returns them as a list
    """
    with open(list_of_ids) as f:
        content = f.readlines()
    content = [x.strip('\n') for x in content]
    return content


def usearch_global(fastq_filename):
    """
    Run vsearch on a fastq file. No limit to # of hits. Minimum 80% identity.
    """
    query_target = fastq_filename
    ref_db = ROOT_PATH + DATABASE_PATH
    p = Popen('vsearch --usearch_global {0} --db {1} '
              ' --id 0.8 --maxaccepts 0 --samout {2}.sam'.format(query_target, ref_db, fastq_filename),
              #wd=os.path.dirname(fastq_filename),
              shell=True,
              executable="/bin/bash")
    p.wait()


def fastq2fasta(fastq_filename):
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


def quality_trim(fastq_filename):
    """
    Quality trimming with BBDuk.
    Quality score threshold of 10 (Brian Bushnell suggestion for Illumina reads).
    Minimum length of 75 to replicate Pal et al. (2016)
    """
    fastq_filename_filtered = fastq_filename.replace('.fastq.gz', '.filtered.fastq.gz')
    filepath_in = fastq_filename
    filepath_out = fastq_filename_filtered

    # Quality-trim to Q10 using Phred algorithm with BBDuk. Trims left and right sides of reads.
    p = Popen('./bbduk.sh -Xmx1g in={0} out={1} qtrim=rl trimq=10 '
              'minlen=75 overwrite=true'.format(filepath_in, filepath_out),
              cwd=BBDUK_PATH,
              shell=True,
              executable='/bin/bash')
    p.wait()


def run(id_list, myfunc):
    """
    Pass in an ID list and apply a given function to all .fastq.gz files
    """
    filepaths_to_process = []
    for dataset_id in id_list:
        dir_contents = os.listdir(MGRAST_DOWNLOAD_PATH + dataset_id)
        for item in dir_contents:
            if item.endswith('.fastq.gz'):
                filepaths_to_process.append(item)

    for filepath in filepaths_to_process:
        myfunc(filepath)


def run_genesippr(id_list):
    """
    TODO: Starts the genesippr docker container then runs against an ID list
    """
    # docker run -it -v /mnt/nas:/mnt/nas genesipprv2
    # python3 genesippr.py /mnt/nas/Forest/MG-RAST_Dataset_Analysis/
        # -s /mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/
        # -t /mnt/nas/Forest/MG-RAST_Dataset_Analysis/db

id_list = retrieve_id_list(ALL_DATASETS)
print(id_list)

# run(id_list, quality_trim)

#fastq2fasta('/mnt/nas/bio_requests/9343/4569599.3/4569599.3.050.upload.filtered.fastq.gz')
usearch_global('/mnt/nas/bio_requests/9343/4569599.3/4569599.3.050.upload.filtered.fasta')