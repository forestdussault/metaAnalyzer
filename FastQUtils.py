import os
import argparse
import time
import glob
import subprocess


class FastQUtils(object):

    def fastq2fasta_bbduk(self):
        """
        Converts .fastq to .fasta. Should be done downstream of quality filtering.
        Should probably deprecate this since the BBTools suite is flexible enough handle compressed/uncompressed filed
        """
        print('\nConverting %s to fasta via --fastq_filter ...\n' % self.fastq_filenames)

        if '.filtered' in self.fastq_filenames:
            fasta_filename = self.fastq_filenames.replace('.fastq.gz', '.fasta')

            # Run VSearch
            p = subprocess.Popen('vsearch '
                                 '--fastq_filter {0} '
                                 '-fastaout {1}'.format(self.fastq_filenames, fasta_filename),
                                 cwd=os.path.dirname(self.fastq_filenames),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()

        else:
            print('No filtered FastQ file provided. Skipping %s' % self.fastq_filenames)
            pass

    @staticmethod
    def run_dedupe(read1, read2):
        """
        Run this after filtering
        http://seqanswers.com/forums/archive/index.php/t-39270.html
        :return:
        """
        # Setup filename
        output = read1.replace('.filtered', '.dereplicated.filtered')
        
        # Run dedupe
        p = subprocess.Popen('dedupe.sh '
                             'in1={0} '
                             'in2={1} '
                             'out={2} '
                             'maxsubs=0 '
                             'ac=f'.format(read1, read2, output),
                             shell=True,
                             executable='/bin/bash')
        p.wait()
        
        # Return the output filename
        return output

    def quality_trim_bbduk(self):
        """
        Quality trimming with BBDuk
        Quality score threshold of 10 (Brian Bushnell suggestion for Illumina reads)
        """
        # Single read handling
        if self.num_reads == 1:
            print('\nRunning BBDuk on %s...\n' % self.fastq_filenames[0])

            # Filename setup
            fastq_filename_filtered = self.fastq_filenames[0].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in = self.fastq_filenames[0]
            filepath_out = fastq_filename_filtered
            stats_out = filepath_out.replace('.fastq.gz', '.paired.stats.txt')

            # Run BBDuk
            p = subprocess.Popen('bbduk.sh '
                                 '-Xmx1g '
                                 'in={0} '
                                 'out={1} '
                                 'qtrim=r '
                                 'trimq=20 '  # Quality-trim to Q20 using Phred algorithm with BBDuk.
                                 'ktrim=r '  # adapter sequences sourced from /bbduk/resources/ should be trimmed from right
                                 'ref={2}/resources/adapters.fa '  # remove adapter sequences
                                 'ktrim=r ' 
                                 'minlen=100 '
                                 'overwrite=true '
                                 'k=23 '  # 23-mers for matching in main portion of read
                                 'tbo '  # BBDuk will internally use BBMerge to trim adapters based on read insert size
                                 'stats={3} '
                                 'threads=8'.format(filepath_in, filepath_out, self.bbduk_dir, stats_out),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()
        
        # Paired read handling
        elif self.num_reads == 2:
            print('\nRunning BBDuk on the provided pair of FastQ files: '
                  '\n%s AND %s' % (self.fastq_filenames[0], self.fastq_filenames[1]))

            # Setup filenames
            fastq_filename_filtered_r1 = self.fastq_filenames[0].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in_r1 = self.fastq_filenames[0]
            filepath_out_r1 = fastq_filename_filtered_r1

            fastq_filename_filtered_r2 = self.fastq_filenames[1].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in_r2 = self.fastq_filenames[1]
            filepath_out_r2 = fastq_filename_filtered_r2

            stats_out = filepath_out_r1.replace('.fastq.gz', '.paired.stats.txt')
            
            # Run BBDuk
            p = subprocess.Popen('bbduk.sh -Xmx1g '
                                 'in1={0} in2={1} '
                                 'out1={2} out2={3} '
                                 'qtrim=r '
                                 'trimq=20 '
                                 'ktrim=r '  # adapter sequences sourced from /bbduk/resources/ should be trimmed from right
                                 'tbo ' 
                                 'tpe '  # trim both read1 and read2
                                 'k=23 '
                                 'ref={4}/resources/adapters.fa '
                                 'minlen=100 '
                                 'overwrite=true '
                                 'stats={5} '
                                 'threads=8'.format(filepath_in_r1, filepath_in_r2, filepath_out_r1, filepath_out_r2,
                                                    self.bbduk_dir, stats_out),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()
            
            # Return the output files
            return filepath_out_r1, filepath_out_r2

    def run_fastqc(self, read1=None, read2=None):
        """
        Run FastQC on paired or single FastQ files
        """
        if self.num_reads == 1:
            # Nasty roundabout way to pass self attributes as default function parameters
            if read1 is None:
                read1 = self.fastq_filenames[0]

            print('\nRunning FastQC on %s...\n' % read1)
            
            # Run FastQC
            p = subprocess.Popen('fastqc {0} --extract'.format(read1),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()
        elif self.num_reads == 2:
            if read1 is None and read2 is None:
                read1 = self.fastq_filenames[0]
                read2 = self.fastq_filenames[1]

            read_list = [read1, read2]

            print('\nRunning FastQC on the provided pair of FastQ files:'
                  '\n%s AND %s...' % (read1, read2))
            
            # Process both reads
            for read in read_list:
                # Run FastQC
                p = subprocess.Popen('fastqc {0} --extract'.format(read),
                                     shell=True,
                                     executable='/bin/bash')
                p.wait()

        # Delete the .zip files produced by FastQC since they have already been extracted into folders at this point
        to_delete = glob.glob(os.path.dirname(self.fastq_filenames[0]) + '/' + '*fastqc.zip')
        for item in to_delete:
            os.remove(item)

    def __init__(self, args):
        self.args = args
        self.fastq_filenames = args.fastq_filenames
        self.fastqc = args.fastqc
        self.fasta = args.fasta
        self.qualitytrim = args.qualitytrim
        self.dereplicate = args.dereplicate
        self.num_reads = len(self.fastq_filenames)

        # Figure out where bbduk is so that we can use the adapter file.
        cmd = 'which bbduk.sh'
        self.bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
        self.bbduk_dir = self.bbduk_dir.split('/')[:-1]
        self.bbduk_dir = '/'.join(self.bbduk_dir)

        if self.num_reads > 2:
            print('\nInvalid number of arguments for .fastq files passed. '
                  'Exiting.')
            quit()

        if self.dereplicate and self.qualitytrim is False:
            print('\nCannot run dereplication unless quality trimming with the -qt flag is specified as well. '
                  'Exiting.')
            quit()

        # Run BBDuk and then FastQC on the fastq pair
        # -qt -dr -fc
        if self.qualitytrim and self.fastqc and self.dereplicate:
            print('\033[92m' + '\033[1m' + '\nRunning BBDuk ==> Dereplication ==> FastQC pipeline... ' + '\033[0m')
            read1_filtered, read2_filtered = self.quality_trim_bbduk()
            interleaved_reads_dedupe = self.run_dedupe(read1_filtered, read2_filtered)
            self.run_fastqc(read1=interleaved_reads_dedupe)
        # -qt -fc
        elif self.qualitytrim and self.fastqc:
            print('\033[92m' + '\033[1m' + '\nRunning BBDuk ==> FastQC pipeline...' + '\033[0m')
            read1_filtered, read2_filtered = self.quality_trim_bbduk()
            self.run_fastqc(read1=read1_filtered, read2=read2_filtered)
        # -qt
        elif self.qualitytrim:
            self.qualitytrim()
        # -fc
        elif self.fastqc:
            self.run_fastqc()

        # -fa
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
    parser.add_argument('-dr', '--dereplicate',
                        default=False,
                        action='store_true',
                        help='Perform dereplication of FastQ with dedupe (BBMap). Important for metagenomic samples.'
                             'Must be run in addition to -qt flag, as this is intended '
                             'to be run post-quality trimming/filtering. ')
    arguments = parser.parse_args()

    x = FastQUtils(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    # Bold green time courtesy of Adam and Andrew
    print('\033[92m' + '\033[1m' + '\nFinished FastQUtils functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')
