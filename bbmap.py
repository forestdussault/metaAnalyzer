import os
import glob
from subprocess import Popen
from config import AMR_DB_PATH, CARB_DB_PATH, BBMAP_TOOLS_PATH


def bbmap_search(fastq_filename, ref_db_fasta):

    if fastq_filename.endswith('.fastq.gz'):
        sam_filename = fastq_filename.replace('.fastq.gz','.sam')
        print('\nWorking with .fastq.gz...')
    elif fastq_filename.endswith('.fastq'):
        sam_filename = fastq_filename.replace('.fastq', '.sam')
        print('\nWorking with .fastq...')
    else:
        print('Invalid input file')
        return None

    base_fastq = os.path.basename(fastq_filename)
    rast_id = base_fastq[:9]
    wd = os.path.dirname(fastq_filename)
    covstats = '{0}/{1}_covstats.txt'.format(wd,rast_id)
    covhist = '{0}/{1}_covhist.txt'.format(wd,rast_id)
    basecov = '{0}/{1}_basecov.txt'.format(wd,rast_id)
    bincov = '{0}/{1}_bincov.txt'.format(wd,rast_id)

    p = Popen('./bbmap.sh in={0} '
              'outm={1} ' # outm to only show mapped alignments
              'ref={2} '
              'nodisk '
              'showprogress=2000000 ' # add progress dotter
              'slow ' # high sensitivity
              'k=12 '
              'trd ' # force old-style cigar strings for compatibility reasons
              'sam=1.3 '
              'covstats={3} ' # produce statistics
              'covhist={4} '
              'basecov={5} '
              'bincov={6}'.format(fastq_filename,sam_filename, ref_db_fasta, covstats, covhist, basecov, bincov),
              cwd=BBMAP_TOOLS_PATH,
              shell=True,
              executable="/bin/bash")
    p.wait()


def metagenome_paths(parent_folder):
    metagenome_folder_list = []
    for filename in glob.glob(parent_folder + '/*/*.fastq.gz', recursive=True):
        metagenome_folder_list.append(filename)
    return metagenome_folder_list


# Grab all of the dataset IDs
metagenome_path_list = metagenome_paths('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes')

# Search all datasets (see bbmap_results_analysis.ipynb for results) with the AMR db
for metagenome_path in metagenome_path_list:
    bbmap_search(metagenome_path,
                 '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/databases/amr/amr_db_original.fasta')