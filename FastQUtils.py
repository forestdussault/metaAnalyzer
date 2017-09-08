import os
import argparse
import time
import glob
from subprocess import Popen

from config import BBMAP_TOOLS_PATH, \
    FASTQC_PATH

class FastQUtils(object):

    def fastq2fasta_bbduk(self):
        '''
        Converts .fastq to .fasta. Should be done downstream of quality filtering.
        '''
        print('\nConverting %s to fasta via --fastq_filter ...\n' % self.fastq_filenames)

        if '.filtered' in self.fastq_filenames:
            fasta_filename = self.fastq_filenames.replace('.fastq.gz', '.fasta')

            p = Popen('vsearch --fastq_filter {0} -fastaout {1}'.format(self.fastq_filenames, fasta_filename),
                      cwd=os.path.dirname(self.fastq_filenames),
                      shell=True,
                      executable='/bin/bash')
            p.wait()

        else:
            print('No filtered fastq file provided. Skipping %s' % self.fastq_filenames)
            pass


    def quality_trim_bbduk(self):
        '''
        Quality trimming with BBDuk
        Quality score threshold of 10 (Brian Bushnell suggestion for Illumina reads)
        Minimum length of 75 to replicate Pal et al. (2016) though all reads should be ~100bp from MiSeq runs
        '''
        if self.num_reads == 1:
            print('\nRunning BBDuk on %s...\n' % self.fastq_filenames[0])

            fastq_filename_filtered = self.fastq_filenames[0].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in = self.fastq_filenames[0]
            filepath_out = fastq_filename_filtered

            # Quality-trim to Q10 using Phred algorithm with BBDuk. Trims left and right sides of reads.
            # Also remove adapter sequences (i.e. Nextera) with the provided /bbmap/resources/adapters.fa file.
            p = Popen('./bbduk.sh -Xmx1g in={0} out={1} qtrim=rl trimq=10 ref=adapters.fa '
                      'minlen=75 overwrite=true'.format(filepath_in, filepath_out),
                      cwd=BBMAP_TOOLS_PATH,
                      shell=True,
                      executable='/bin/bash')
            p.wait()

        elif self.num_reads == 2:
            print('\nRunning BBDuk on the provided pair of FastQ files: '
                  '\n%s AND %s' % (self.fastq_filenames[0],self.fastq_filenames[1]))

            fastq_filename_filtered_r1 = self.fastq_filenames[0].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in_r1 = self.fastq_filenames[0]
            filepath_out_r1 = fastq_filename_filtered_r1

            fastq_filename_filtered_r2 = self.fastq_filenames[1].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in_r2 = self.fastq_filenames[1]
            filepath_out_r2 = fastq_filename_filtered_r2

            # Pairs should always be kept together with BBDuk, hence usage of in1, in2, out1, out2
            p = Popen('./bbduk.sh -Xmx1g '
                      'in1={0} in2={1} '
                      'out1={2} out2={3} '
                      'qtrim=rl trimq=10 '
                      'tbo ' 
                      'tpe '
                      'ref=./resources/adapters.fa '
                      'minlen=75 overwrite=true'.format(filepath_in_r1, filepath_in_r2, filepath_out_r1, filepath_out_r2),
                      cwd=BBMAP_TOOLS_PATH,
                      shell=True,
                      executable='/bin/bash')
            p.wait()
            return filepath_out_r1, filepath_out_r2

    def run_fastqc(self, read1=None, read2=None):
        if self.num_reads == 1:
            if read1 == None:
                read1 = self.fastq_filenames[0]

            print('\nRunning FastQC on %s...\n' % read1)
            p = Popen('./fastqc {0} --extract'.format(read1),
                      cwd=FASTQC_PATH,
                      shell=True,
                      executable='/bin/bash')
            p.wait()
        elif self.num_reads == 2:
            if read1 == None and read2 == None:
                read1 = self.fastq_filenames[0]
                read2 = self.fastq_filenames[1]

            read_list = [read1,read2]

            print('\nRunning FastQC on the provided pair of FastQ files:'
                  '\n%s AND %s...' % (read1,read2))
            for read in read_list:
                p = Popen('./fastqc {0} --extract'.format(read),
                          cwd=FASTQC_PATH,
                          shell=True,
                          executable='/bin/bash')
                p.wait()

        # Delete the .zip files produced by fastqc since they have already been extracted into the folder at this point
        to_delete = glob.glob(os.path.dirname(self.fastq_filenames[0]) + '/' + '*fastqc.zip')
        for item in to_delete:
            print(item)
            os.remove(item)

    def __init__(self, args):
        self.args = args
        self.fastq_filenames = args.fastq_filenames
        self.fastqc = args.fastqc
        self.fasta = args.fasta
        self.qualitytrim = args.qualitytrim
        self.num_reads = len(self.fastq_filenames)

        if self.num_reads > 2:
            print('\nInvalid number of arguments for .fastq files passed. Exiting.')
            quit()

        if self.qualitytrim and self.fastqc:
            print('\nRunning BBDuk and then FastQC afterwards...')
            read1_filtered, read2_filtered = self.quality_trim_bbduk()
            self.run_fastqc(read1=read1_filtered, read2=read2_filtered)
        elif self.qualitytrim:
            self.qualitytrim()
        elif self.fastqc:
            self.run_fastqc()

        if self.fasta:
            self.fasta()

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_filenames',
                           help='Path to the FastQ file you would like to perform QC functions on. '
                                'Enter two reads for pairs i.e. /path/read1.fastq /path/read2.fastq.', nargs='*')
    parser.add_argument('-fa', '--fasta',
                        default=False,
                        action='store_true',
                        help='Convert .fastq to .fasta format.')
    parser.add_argument('-qt', '--qualitytrim',
                        default=False,
                        action='store_true',
                        help='Perform trimming with BBDuk on target FastQ file.'
                             'If the parameter for running FastQC is also passed, the trimmed files will be passed to'
                             'FastQC after filtering is complete.')
    parser.add_argument('-fc', '--fastqc',
                        default=False,
                        action='store_true',
                        help='Perform analysis on FastQ file with FastQC. '
                             'If the parameter for running BBDuk quality trimming is also passed, FastQC will be run '
                             'on the filtered files instead of the originals.')
    arguments = parser.parse_args()

    x = FastQUtils(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    # Bold green time courtesy of Adam and Andrew
    print('\033[92m' + '\033[1m' + '\nFinished FastQUtils functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')