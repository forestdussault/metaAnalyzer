import os
import glob
import subprocess


def bbmap_search_sensitive(fastq_filename, ref_db_fasta):

    if fastq_filename.endswith('.fastq.gz'):
        sam_filename = fastq_filename.replace('.fastq.gz', '.sam')
    elif fastq_filename.endswith('.fastq'):
        sam_filename = fastq_filename.replace('.fastq', '.sam')
    else:
        print('Invalid input file')
        return None

    base_fastq = os.path.basename(fastq_filename)
    rast_id = base_fastq[:9]
    wd = os.path.dirname(fastq_filename)
    covstats = '{0}/{1}_covstats.txt'.format(wd, rast_id)
    covhist = '{0}/{1}_covhist.txt'.format(wd, rast_id)
    basecov = '{0}/{1}_basecov.txt'.format(wd, rast_id)
    bincov = '{0}/{1}_bincov.txt'.format(wd, rast_id)

    p = subprocess.Popen('bbmap.sh in={0} '
                         'outm={1} '  # outm to only show mapped alignments
                         'ref={2} '
                         'nodisk '
                         'showprogress=2000000 '  # add progress dotter
                         'slow '  # high sensitivity
                         'k=12 '
                         'trd '  # force old-style cigar strings for compatibility reasons
                         'sam=1.3 '
                         'covstats={3} '  # produce statistics
                         'covhist={4} '
                         'basecov={5} '
                         'bincov={6}'.format(fastq_filename,
                                             sam_filename,
                                             ref_db_fasta,
                                             covstats,
                                             covhist,
                                             basecov,
                                             bincov),
                         shell=True,
                         executable="/bin/bash")
    p.wait()


def bbmap_search_fast(fastq_filename, ref_db_fasta):

    if fastq_filename.endswith('.fastq.gz'):
        sam_filename_aligned = fastq_filename.replace('.fastq.gz', '.aligned.sam')
        sam_filename_unaligned = fastq_filename.replace('.fastq.gz', '.unaligned.sam')
    elif fastq_filename.endswith('.fastq'):
        sam_filename_aligned = fastq_filename.replace('.fastq', '.aligned.sam')
        sam_filename_unaligned = fastq_filename.replace('.fastq', '.unaligned.sam')
    else:
        print('Invalid input file')
        return None

    base_fastq = os.path.basename(fastq_filename)
    rast_id = base_fastq[:9]
    wd = os.path.dirname(fastq_filename)
    covstats = '{0}/{1}_covstats.txt'.format(wd, rast_id)
    covhist = '{0}/{1}_covhist.txt'.format(wd, rast_id)
    basecov = '{0}/{1}_basecov.txt'.format(wd, rast_id)
    bincov = '{0}/{1}_bincov.txt'.format(wd, rast_id)
    scafstats = '{0}/{1}_scafstats.txt'.format(wd, rast_id)

    p = subprocess.Popen('bbmap.sh in={0} '
                         'outm={1} '
                         'outu={2}'  # outm to only show mapped alignments
                         'ref={3} '
                         'nodisk '
                         'showprogress=2000000 '  # add progress dotter
                         'fast ' 
                         'trd '  # force old-style cigar strings for compatibility reasons
                         'sam=1.3 '
                         'covstats={4} '  # produce statistics
                         'covhist={5} '
                         'basecov={6} '
                         'bincov={7} '
                         'scafstats={8} '
                         'bs=bs.sh; sh bs.sh'  # sort the BAM file
                         ''.format(fastq_filename,
                                   sam_filename_aligned,
                                   sam_filename_unaligned,
                                   ref_db_fasta,
                                   covstats,
                                   covhist,
                                   basecov,
                                   bincov,
                                   scafstats),
                         shell=True,
                         executable="/bin/bash")
    p.wait()


def bbmap_search_fast_paired(fastq_r1, fastq_r2, ref_db_fasta):
    """
    Mapping according to B. Bushnell's recommendation for parameters for
    "very high precision and lower sensitivity, as when removing contaminant
    reads specific to a genome without risking false-positives" (see BBMap Guide)
    """
    if fastq_r1.endswith('.fastq.gz') and fastq_r2.endswith('.fastq.gz'):
        sam_output = fastq_r1.replace('.fastq.gz', '.aligned.sam')
    elif fastq_r1.endswith('.fastq'):
        sam_output = fastq_r1.replace('.fastq', '.aligned.sam')
    else:
        print('Invalid input files')
        return None

    base_fastq = os.path.basename(fastq_r1)
    wd = os.path.dirname(fastq_r1)
    covstats = '{0}/{1}_covstats.txt'.format(wd, base_fastq)
    covhist = '{0}/{1}_covhist.txt'.format(wd, base_fastq)
    bincov = '{0}/{1}_bincov.txt'.format(wd, base_fastq)
    scafstats = '{0}/{1}_scafstats.txt'.format(wd, base_fastq)
    basecov = '{0}/{1}_basecov.txt'.format(wd, base_fastq)

    p = subprocess.Popen('bbmap.sh in1={0} '
                         'in2={1} '
                         'out={2} '
                         'ref={3} '
                         'nodisk '  # don't write ref. index to disk -> keep index in memory
                         'showprogress=2000000 '  # add progress dotter
                         'fast '  # all of the following parameters up til covstats are for high precision
                         'minratio=0.9 '
                         'maxindel=3 '
                         'bwr=0.16 '
                         'minhits=2 '
                         'qtrim=r '
                         'trimq=10 '
                         'untrim '
                         'idtag '
                         'printunmappedcount '
                         'kfilter=25 '
                         'maxsites=1 '
                         'k=14 '
                         'covstats={4} '  # produce statistics
                         'covhist={5} '
                         'bincov={6} '
                         'scafstats={7} '
                         'basecov={8} '
                         'bs=bs.sh; sh bs.sh'  # sort the BAM file
                         ''.format(fastq_r1, fastq_r2, sam_output, ref_db_fasta,
                                   covstats, covhist, bincov, scafstats, basecov),
                         shell=True,
                         executable="/bin/bash")
    p.wait()


def metagenome_paths(parent_folder):
    metagenome_folder_list = []
    for filename in glob.glob(parent_folder + '/*/*.fastq.gz', recursive=True):
        metagenome_folder_list.append(filename)
    return metagenome_folder_list


# Grab all of the dataset IDs
# metagenome_path_list = metagenome_paths('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes')

# Search all datasets (see bbmap_results_analysis.ipynb for results) with the AMR db
# for metagenome_path in metagenome_path_list:
#     bbmap_search_sensitive(metagenome_path,
#                  CARB_DB_PATH)

bbmap_search_fast_paired('/mnt/scratch/Forest/SRA_carrot_project/metagenomes/qualimap_test/SRR3747715_1.filtered.fastq.gz',
                         '/mnt/scratch/Forest/SRA_carrot_project/metagenomes/qualimap_test/SRR3747715_2.filtered.fastq.gz',
                         '/mnt/scratch/Forest/SRA_carrot_project/genomes/cow/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna')

# Carrot
# '/mnt/scratch/Forest/SRA_carrot_project/genomes/carrot/GCF_001625215.1_ASM162521v1_genomic.fna')