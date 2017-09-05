import os
import docker
import shutil

from glob import glob
from subprocess import Popen

from config import AMR_DB_PATH, \
    MGRAST_DOWNLOAD_PATH, \
    BBMAP_TOOLS_PATH


def retrieve_id_list(list_of_ids):
    """
    Parses a text file with individual IDs on newlines and returns them as a list
    """
    with open(list_of_ids) as f:
        content = f.readlines()
    content = [x.strip('\n') for x in content]
    return content


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


def quality_trim(fastq_filename):
    """
    Quality trimming with BBDuk
    Quality score threshold of 10 (Brian Bushnell suggestion for Illumina reads)
    Minimum length of 75 to replicate Pal et al. (2016)
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


def create_genesippr_container():
    """
    Starts the genesippr docker container
    """

    # Instantiate client to talk to Docker daemon
    client = docker.from_env()

    container = client.containers.run("genesipprv2:latest",
                                      detach=True,
                                      stdin_open=True,
                                      tty=True,
                                      volumes={'/mnt/nas': '/mnt/nas'})

    return container


def run_genesippr(container, input_directory, sequence_path, target_path):
    """
    Takes a genesippr container and runs the program
    """
    a = container.exec_run('python3 /geneSipprV2/sipprverse/genesippr/genesippr.py {0}' 
                           ' -s {1}' 
                           ' -t {2}' 
                           ' --detailedReports'.format(input_directory, sequence_path, target_path))
    print("{}".format(a.strip().decode('UTF-8')))


def setup_symlinks():
    """
    One off function. Used this to create symlinks in my analysis folder
    """
    subdirs = glob('/mnt/nas/bio_requests/9343/*/')
    dirs_to_create = []
    sym_link_refs = {}
    for dir in subdirs:
        temp_file_list = glob(dir + '*.filtered.fastq.gz')
        if len(temp_file_list) > 0:
            sym_link_refs[temp_file_list[0]] = temp_file_list[0].replace(
                '/bio_requests/9343/',
                '/Forest/MG-RAST_Dataset_Analysis/metagenomes/')
            temp_file_list[0] = temp_file_list[0].replace('/bio_requests/9343/',
                                                          '/Forest/MG-RAST_Dataset_Analysis/metagenomes/')
            dirs_to_create.append(temp_file_list[0])

    # Create folders
    for dir in dirs_to_create:
        if not os.path.exists(dir):
            print('mkdir ' + dir[:-38])
            os.makedirs(dir[:-38])

    # Create symlinks
    for key, value in sym_link_refs.items():
        os.symlink(key, value)

def mass_sipp(id_list):
    # Get list of folders to run genesippr on
    id_list = os.listdir('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes')

    # Start the container
    container = create_genesippr_container()

    # Send commands to the container
    for id in id_list:
        print("\nRunning genesippr on {}...".format(id))
        run_genesippr(container=container,
                      input_directory='/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/{}'.format(id),
                      sequence_path='/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/{}'.format(id),
                      target_path='/mnt/nas/Forest/MG-RAST_Dataset_Analysis/db')

def clean_folder(parent_folder):
    # Delete everything that isn't .filtered.fastq.gz (including folders)
    for item in os.listdir(parent_folder):
        if item.endswith('.filtered.fastq.gz'):
            print('Skipping %s ...' % item)
        else:
            print('Deleting %s ...' % item)

            if os.path.isdir(parent_folder + '/' + item):
                shutil.rmtree(parent_folder + '/' +item)
            else:
                os.remove(parent_folder + '/' +item)
