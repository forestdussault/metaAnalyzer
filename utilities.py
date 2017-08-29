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


# Probably don't need to do this again
# def delete_extraneous_files(mgrast_download_folder):
#     """
#     Deletes everything that is not a .fastq.gz file in the passed in folder
#     """
#     for file in os.listdir(mgrast_download_folder):
#         if 'upload.fastq.gz' in file:
#             print('Keeping %s ...' % file)
#             pass
#         else:
#             print('Deleting %s ...' % file)
#
#             try:
#                 os.remove((mgrast_download_folder + '/' + file))
#             except:
#                 print('Could not delete %s. Skipping.' % file)
#                 pass


def usearch_global(mgrast_id):
    """
    Run vsearch on a fastq file. No limit to # of hits. Minimum 90% identity.
    """
    p = Popen('vsearch --usearch_global %s --db %s '
              ' --id 0.9 --maxaccepts 0 --samout samout.txt' % (ROOT_PATH + mgrast_id, ROOT_PATH + DATABASE_PATH),
              #wd=self.mg_rast_tools,
              shell=True,
              executable="/bin/bash")
    p.wait()


def fastq_filter(mgrast_id):
    """
    Converts .fastq to .fasta
    """
    print('\nConverting %s to fasta via --fastq_filter ...' % mgrast_id)

    if '.filtered' in mgrast_id:
        mgrast_id_fasta = mgrast_id.replace('.fastq.gz','.fasta')

        try:
            p = Popen('vsearch --fastq_filter {0} -fastaout {1}'.format(mgrast_id,mgrast_id_fasta),
                cwd = MGRAST_DOWNLOAD_PATH+mgrast_id[:9],
                shell = True,
                executable = '/bin/bash')
            p.wait()
        except:
            print('Could not start fastq -> fasta conversion.')
    else:
        pass


def quality_trim(fastq_filename):
    """
    Does some quality trimming with BBDuk
    """
    fastq_filename_filtered = fastq_filename.replace('.fastq.gz','.filtered.fastq.gz')
    mgrast_id = fastq_filename[:9]
    filepath_in = MGRAST_DOWNLOAD_PATH + mgrast_id + '/' + fastq_filename
    filepath_out = MGRAST_DOWNLOAD_PATH + mgrast_id + '/' + fastq_filename_filtered
    try:
        #Quality-trim to Q10 using Phred algorithm with BBDuk. Trims left and right sides of reads.
        p = Popen('./bbduk.sh -Xmx1g in={0} out={1} qtrim=rl trimq=10 minlen=75 overwrite=true'.format(filepath_in, filepath_out),
            cwd = BBDUK_PATH,
            shell = True,
            executable = '/bin/bash')
        p.wait()
    except:
        print('Could not start BBDuk quality trimming.')


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

id_list = retrieve_id_list(ALL_DATASETS)
print(id_list)

run(id_list, quality_trim)



