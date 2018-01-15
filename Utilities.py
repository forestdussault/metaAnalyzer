import os
import shutil
from glob import glob
from config import MGRAST_DOWNLOAD_PATH


def retrieve_id_list(list_of_ids):
    """
    Parses a text file with individual IDs on newlines and returns them as a list
    """
    with open(list_of_ids) as f:
        content = f.readlines()
    content = [x.strip('\n') for x in content]
    return content


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


def setup_symlinks():
    """
    One off function. Used this to create symlinks in my analysis folder
    """
    subdirs = glob('/mnt/nas/bio_requests/9343/*/')
    dirs_to_create = []
    sym_link_refs = {}
    for dir in subdirs:
        temp_file_list = glob(dir + '*.upload.filtered.fastq.gz')

        for item in temp_file_list:
            if 'dereplicated' in item:
                temp_file_list.remove(item)

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
            try:
                os.makedirs(dir[:-38])
            except:
                pass

    # Create symlinks
    for key, value in sym_link_refs.items():
        os.symlink(key, value)


def clean_folder(parent_folder):
    # Delete everything that isn't .fastq.gz (including folders)
    for item in os.listdir(parent_folder):
        if item.endswith('.fastq.gz'):
            print('Skipping %s ...' % item)
        else:
            print('Deleting %s ...' % item)

            if os.path.isdir(parent_folder + '/' + item):
                shutil.rmtree(parent_folder + '/' + item)
            else:
                os.remove(parent_folder + '/' + item)


def move_results(parent_folder, name):
    # Move files that aren't .fastq.gz or folders into a new named folder
    try:
        print('Creating directory ' + parent_folder + '/' + name)
        os.mkdir(parent_folder + '/' + name)
    except FileExistsError:
        return None
    for item in os.listdir(parent_folder):
        if item.endswith('.filtered.fastq.gz') or os.path.isdir(parent_folder + '/' + item):
            print('Skipping %s ...' % item)
        else:
            print('Moving %s ...' % item)
            shutil.move((parent_folder + '/' + item),(parent_folder + '/' + name + '/' + item))

# for item in os.listdir('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes'):
#     move_results('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/' + item, 'bbmap_carbapenemase_amr_db_results')
#

# setup_symlinks()