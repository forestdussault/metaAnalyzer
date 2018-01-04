"""
1. Delete extraneous files from  /mnt/nas/bio_requests/9343/{dataset_id}

2. Perform filtering with bbduk on contents of download folder to create .filtered files

3. Create symlinks from bio_request folder of downloads to the analysis folder (/mnt/nas/Forest/MG-RAST_Dataset_Analysis)

"""

import os
import glob
import subprocess


def retrieve_id_list(idfile):
    with open(idfile) as f:
        content = f.readlines()
    content = [x.strip('\n') for x in content]
    return content


# Get IDs to work with
id_list = retrieve_id_list('/mnt/nas/bio_requests/9343/datasets.txt')

# Get subdirectories to crawl through
subdirs = glob.glob('/mnt/nas/bio_requests/9343/*/')


# Filter all raw sequences with FastQUtils
for directory in subdirs:
    print('\nWorking on:\t{}...'.format(directory))
    file_list = os.listdir(directory)
    fastq = None
    for file in file_list:
        if file.endswith('.upload.fastq.gz'):
            fastq = os.path.join(os.path.dirname(directory), file)

    if fastq is None:
        print('Skipping:\t{}'.format(fastq))
    else:
        cmd = 'python3 /home/dussaultf/PycharmProjects/mgrast_analysis/FastQUtils.py {fastq_filepath} ' \
              '-qt -dr '.format(fastq_filepath=fastq)

        p = subprocess.Popen(cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()
        print('Done with:\t{}'.format(fastq))

# # Delete extraneous files
# for directory in subdirs:
#     for item in os.listdir(directory):
#         # Only target files that end with .fastq.gz
#         if item.endswith('.fastq.gz'):
#             # Ignore raw source files (.upload.) and filtered files (.filtered.)
#             if item.endswith('.upload.fastq.gz'):
#                 print('Keeping {item}'.format(item=item))
#             else:
#                 print('Deleting {item}'.format(item=item))
#                 os.remove(os.path.join(directory, item))
#         else:

#             pass