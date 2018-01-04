import subprocess
import os
import glob


def bbmap_sensitive(fastq, refdb, db_prepend):

    if fastq.endswith('.fastq.gz'):
        sam_filename = fastq.replace('.fastq.gz', '_{}.sam'.format(db_prepend))
    elif fastq.endswith('.fastq'):
        sam_filename = fastq.replace('.fastq', '_{}.sam'.format(db_prepend))
    else:
        print('Invalid input file')
        return None

    rast_id = os.path.basename(fastq)[:9]

    workdir = os.path.dirname(fastq)

    try:
        os.makedirs(os.path.join(workdir, db_prepend))
    except OSError:
        print('Skipping {}'.format(fastq))
        return None

    workdir = os.path.join(workdir, db_prepend)
    sam_filename = os.path.basename(sam_filename)
    sam_filename = os.path.join(workdir, sam_filename)


    covstats = '{0}/{1}_{2}_covstats.txt'.format(workdir, rast_id, db_prepend)
    covhist = '{0}/{1}_{2}_covhist.txt'.format(workdir, rast_id, db_prepend)
    basecov = '{0}/{1}_{2}_basecov.txt'.format(workdir, rast_id, db_prepend)
    bincov = '{0}/{1}_{2}_bincov.txt'.format(workdir, rast_id, db_prepend)

    p = subprocess.Popen('bbmap.sh in={0} '
                         'outm={1} '  # outm to only show mapped alignments
                         'ref={2} '
                         'nodisk '
                         'slow '  # high sensitivity
                         'trd '  # force old-style cigar strings for compatibility reasons
                         'sam=1.3 '
                         'covstats={3} '  # produce statistics
                         'covhist={4} '
                         'basecov={5} '
                         'bincov={6}'.format(fastq,
                                             sam_filename,
                                             refdb,
                                             covstats,
                                             covhist,
                                             basecov,
                                             bincov),
                         shell=True,
                         executable="/bin/bash")
    p.wait()

# amr_db = '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/databases/amr/amr_db_original.fasta'
# carb_db = '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/databases/carbapenemase/carbapenemase_fasta_files/' \
#           'carbapenemase_db/carbapenemase_merged.fasta'
# mem_db = '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/databases/amr/MEMA1_MEMB2.fasta'
sul4_db = '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/databases/sul4/sul4.fasta'

wd = '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes'
work_dirs = os.listdir(wd)

for item in work_dirs:
    path = os.path.join(wd, item)
    fastq = glob.glob(os.path.join(path, '*.filtered.fastq.gz'))
    for item in fastq:
        if 'dereplicated' in item:
            fastq.remove(item)
    bbmap_sensitive(fastq[0], sul4_db, 'sul4_db')
    print('Done with {}\n'.format(os.path.basename(fastq[0])))
