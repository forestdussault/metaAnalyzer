import os
import subprocess
import time
import argparse


class BBMapper(object):

    def bbmap_sensitive(self):

        if self.read1.endswith('.fastq.gz'):
            sam_filename = self.read1.replace('.fastq.gz', '.sam')
        elif self.read1.endswith('.fastq'):
            sam_filename = self.read1.replace('.fastq', '.sam')
        else:
            print('Invalid input file')
            return None

        rast_id = os.path.basename(self.read1)[:9]

        covstats = '{0}/{1}_covstats.txt'.format(self.workdir, rast_id)
        covhist = '{0}/{1}_covhist.txt'.format(self.workdir, rast_id)
        basecov = '{0}/{1}_basecov.txt'.format(self.workdir, rast_id)
        bincov = '{0}/{1}_bincov.txt'.format(self.workdir, rast_id)

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
                             'bincov={6}'.format(self.read1,
                                                 sam_filename,
                                                 self.refdb,
                                                 covstats,
                                                 covhist,
                                                 basecov,
                                                 bincov),
                             shell=True,
                             executable="/bin/bash")
        p.wait()

    def bbmap_fast_paired(self):
        """
        Mapping according to B. Bushnell's recommendation for parameters for
        "very high precision and lower sensitivity, as when removing contaminant
        reads specific to a genome without risking false-positives" (see BBMap Guide)
        """
        if self.read1.endswith('.fastq.gz') and self.read2.endswith('.fastq.gz'):
            sam_output = self.read1.replace('.fastq.gz', '.aligned.sam')
        elif self.read1.endswith('.fastq'):
            sam_output = self.read1.replace('.fastq', '.aligned.sam')
        else:
            print('Invalid input files')
            return None

        base_fastq = os.path.basename(self.read1)

        covstats = '{0}/{1}_covstats.txt'.format(self.workdir, base_fastq)
        covhist = '{0}/{1}_covhist.txt'.format(self.workdir, base_fastq)
        bincov = '{0}/{1}_bincov.txt'.format(self.workdir, base_fastq)
        scafstats = '{0}/{1}_scafstats.txt'.format(self.workdir, base_fastq)
        basecov = '{0}/{1}_basecov.txt'.format(self.workdir, base_fastq)

        p = subprocess.Popen('bbmap.sh '
                             '-Xmx26g '  # pass to Java for memory detection.
                             'in1={0} '
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
                             ''.format(self.read1, self.read2, sam_output, self.refdb,
                                       covstats, covhist, bincov, scafstats, basecov),
                             shell=True,
                             executable="/bin/bash")

        p.wait()

    def run_seal(self, pergene=None):
        """
        Generate abundance measurements
        """
        print('\nRunning Seal...')
        refnames = '' #  PEP8 compliance nonsense
        if pergene is None:
            refnames = 't'
        elif pergene is True:
            refnames = 'f'

        base_fastq = os.path.basename(self.read1)

        stats = '{0}/{1}_sealstats.txt'.format(self.workdir, base_fastq)
        rpkm = '{0}/{1}_sealrpkm.txt'.format(self.workdir, base_fastq)

        p = subprocess.Popen('seal.sh '
                             'in1={0} '
                             'in2={1} '
                             'ref={2} '
                             'stats={3} '
                             'rpkm={4} '
                             'ambig=random '
                             'overwrite=true '
                             'refnames={5}'  # per ref. file instead of per target in ref. (set this to true for host DNA %)
                             ''.format(self.read1, self.read2, self.refdb, stats, rpkm, refnames),
                             shell=True,
                             executable="/bin/bash")
        p.wait()

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nBBMAPPER' + '\033[0m')

        # Get args
        self.args = args

        # Required input
        self.fastq_filenames = args.fastq_filenames
        self.refdb = args.refdb

        # Flags
        self.pergene = args.pergene
        self.seal = args.seal
        self.pairedsensitive = args.pairedsensitive

        # Metadata
        self.num_reads = len(self.fastq_filenames)
        self.workdir = os.path.dirname(self.fastq_filenames[0])

        # Read setup
        if self.num_reads == 1:
            self.read1 = self.fastq_filenames[0]
        elif self.num_reads == 2:
            self.read1 = self.fastq_filenames[0]
            self.read2 = self.fastq_filenames[1]
        else:
            print('Invalid number of reads entered. Exiting.')
            quit()

        # BBDuk directory grab
        cmd = 'which bbduk.sh'
        self.bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
        self.bbduk_dir = self.bbduk_dir.split('/')[:-1]
        self.bbduk_dir = '/'.join(self.bbduk_dir)

        # User input validation
        if self.pergene and self.seal:
            print('Specified too many analysis types. Quitting.')
            quit()
        elif self.pergene and self.seal is False:
            print('Must specify -s flag in order to use -pg flag. Quitting.')
            quit()

        # Run methods
        if self.seal and self.pergene:
            self.run_seal(pergene=self.pergene)
        elif self.seal:
            self.run_seal()
        elif self.pairedsensitive:
            self.bbmap_fast_paired()

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_filenames',
                        help='Path to the FastQ file(s) you would like to perform BBMap functions on. '
                             'Enter two reads for pairs i.e. /path/read1.fastq /path/read2.fastq.',
                        nargs='*')
    parser.add_argument('-db', '--refdb',
                        help='Specify the path to the reference database you would like to query against.',
                        required=True)
    parser.add_argument('-s', '--seal',
                        default=False,
                        action='store_true',
                        help='Runs seal.sh against the input read(s). '
                             'Provides abundance estimates for gene targets or the reference DB treated as a whole. '
                             'By default this will provide only one abundance value'
                             'i.e. the abundance in % of the whole DB in your sample '
                             'Add the -pg flag to provide abundance estimates per sequence in your reference DB.')
    parser.add_argument('-pg', '--pergene',
                        default=False,
                        action='store_true',
                        help='Modifier for the --seal flag '
                             'Use this if your reference DB is a small set of genes of interest.')
    parser.add_argument('-ps', '--pairedsensitive',
                        default=False,
                        action='store_true',
                        help='Runs bbmap.sh against the input read(s). '
                             'Very high precision and lower sensitivity BBMap search against reference DB.'
                             'For removing contaminant reads specific to a metagenome.')

    arguments = parser.parse_args()

    x = BBMapper(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    # Bold green time courtesy of Adam and Andrew
    print('\033[92m' + '\033[1m' + '\nFinished BBMapper functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')


# AMR DB
# '/mnt/nas/Forest/MG-RAST_Dataset_Analysis/databases/amr/amr_db_original.fasta'

# Carrot
# '/mnt/scratch/Forest/SRA_carrot_project/genomes/carrot/GCF_001625215.1_ASM162521v1_genomic.fna'

# Cow
# '/mnt/scratch/Forest/SRA_carrot_project/genomes/cow/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna'
