import os
import argparse
import time
import glob
import subprocess
import shutil


def run_merge(read1, read2):
    print('\nMerging {} and {}...'.format(read1, read2))

    # Strip pair number
    output_filename = read1.replace('_1', '')

    # Set final filename
    output_filename = output_filename.replace('.fastq.gz', '.merged.fastq.gz')

    p = subprocess.Popen('bbmerge.sh '
                         'in1={} '
                         'in2={} '
                         'out={} '
                         'overwrite=true '
                         ''.format(read1, read2, output_filename),
                         shell=True,
                         executable="/bin/bash")
    p.wait()
    return output_filename


class FastQUtils(object):
    def quality_trim_bbduk(self):
        """
        Quality trimming with BBDuk
        """
        # Single read handling
        if self.num_reads == 1:
            print('\nRunning BBDuk on %s...\n' % self.fastq_filenames[0])

            # Filename setup
            fastq_filename_filtered = self.fastq_filenames[0].replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in = self.fastq_filenames[0]
            filepath_out = fastq_filename_filtered
            stats_out = filepath_out.replace('.fastq.gz', '.reference.stats.txt')

            # Run BBDuk
            p = subprocess.Popen('bbduk.sh '
                                 '-Xmx1g '
                                 'in={0} '
                                 'out={1} '
                                 'qtrim=r '
                                 'trimq=20 '  # Quality-trim to Q20 using Phred algorithm with BBDuk.
                                 'ktrim=r '  # adapter sequences sourced from bbduk should be trimmed from right
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

        elif self.num_reads == 2:
            # Paired read handling
            read1 = self.fastq_filenames[0]
            read2 = self.fastq_filenames[1]

            print('\nRunning BBDuk on the provided pair of FastQ files: '
                  '\n%s AND %s' % (read1, read2))

            # Setup filenames
            fastq_filename_filtered_r1 = read1.replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in_r1 = read1
            filepath_out_r1 = fastq_filename_filtered_r1

            fastq_filename_filtered_r2 = read2.replace('.fastq.gz', '.filtered.fastq.gz')
            filepath_in_r2 = read2
            filepath_out_r2 = fastq_filename_filtered_r2

            stats_out = filepath_out_r1.replace('.fastq.gz', '.paired.stats.txt')

            # Run BBDuk
            p = subprocess.Popen('bbduk.sh -Xmx4g '
                                 'in1={0} in2={1} '
                                 'out1={2} out2={3} '
                                 'qtrim=r '
                                 'trimq=20 '
                                 'ktrim=r '  # adapter sequences sourced from bbduk should be trimmed from right
                                 'tbo '
                                 'tpe '  # trim both read1 and read2
                                 'k=23 '
                                 'ref={4}/resources/adapters.fa '
                                 'minlen=100 '
                                 'overwrite=true '
                                 'stats={5}'
                                 ''.format(filepath_in_r1, filepath_in_r2, filepath_out_r1, filepath_out_r2,
                                           self.bbduk_dir, stats_out),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()

            # Return the output files
            return filepath_out_r1, filepath_out_r2

    def run_dedupe(self, read1=None, read2=None):
        """
        Run this after filtering. Deduplicates sequences. Requires lots of memory.
        http://seqanswers.com/forums/archive/index.php/t-39270.html
        """

        # Single read
        if self.num_reads == 1:
            print('\nRunning dedupe on {0}...\n'.format(read1))

            # Setup output
            dedupe_out_1 = read1.replace('.filtered', '.dereplicated.filtered')

            # Run dedupe
            p = subprocess.Popen('dedupe.sh '
                                 'in={0} '
                                 'out={1} '
                                 'maxsubs=0 '
                                 'k=31 '
                                 'ac=f '
                                 'overwrite=true'.format(read1, dedupe_out_1),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()

            return dedupe_out_1

        elif self.num_reads == 2:

            print('\nRunning dedupe on {0} and {1}...\n'.format(read1, read2))

            # Setup temp output
            interleaved_out = read1.replace('.filtered', '.temp.interleaved')

            # Run dedupe
            p = subprocess.Popen('dedupe.sh '
                                 '-Xmx16g '
                                 'in1={0} '
                                 'in2={1} '
                                 'out={2} '
                                 'maxsubs=0 '
                                 'k=31 '
                                 'ac=f '
                                 'overwrite=true'.format(read1, read2, interleaved_out),
                                 shell=True,
                                 executable='/bin/bash')
            p.wait()

            # Setup final output
            dedupe_out_1 = read1.replace('.filtered', '.dereplicated.filtered')
            dedupe_out_2 = read2.replace('.filtered', '.dereplicated.filtered')

            # Deinterleave the output file
            print('\nDeinterleaving...')
            p2 = subprocess.Popen('reformat.sh '
                                  'in={0} '
                                  'out1={1} '
                                  'out2={2} '
                                  'overwrite=true'.format(interleaved_out, dedupe_out_1, dedupe_out_2),
                                  shell=True,
                                  executable='/bin/bash')
            p2.wait()

            # Cleanup
            os.remove(interleaved_out)

            # Return the output filenames
            return dedupe_out_1, dedupe_out_2

    def run_fastqc(self, read1=None, read2=None):
        """
        Run FastQC on paired or single FastQ files
        """
        unfiltered = False
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
                unfiltered = True
                read1 = self.fastq_filenames[0]
                read2 = self.fastq_filenames[1]

            read_list = [read1, read2]

            print('\nRunning FastQC on the provided pair of FastQ files:'
                  '\n%s and %s...' % (read1, read2))
            
            # Process both reads
            for read in read_list:
                # Run FastQC
                p = subprocess.Popen('fastqc {0} '
                                     '--extract'.format(read),
                                     shell=True,
                                     executable='/bin/bash')
                p.wait()

        # Delete the .zip files produced by FastQC since they have already been extracted at this point
        to_delete = glob.glob(self.workdir + '/' + '*fastqc.zip')
        for item in to_delete:
            os.remove(item)

        # Create a directory for the FastQC results
        if unfiltered:
            fastqc_foldername = '/FastQC_unfiltered'
        else:
            fastqc_foldername = '/FastQC_filtered'

        try:
            os.mkdir(self.workdir + fastqc_foldername)
        except OSError:
            shutil.rmtree(self.workdir + fastqc_foldername)
            os.mkdir(self.workdir + fastqc_foldername)

        # Move all _fastqc files/folders to the FastQC directory
        to_move = glob.glob(self.workdir + '/' + '*_fastqc*')

        for file in to_move:
            shutil.move((self.workdir + '/' + os.path.basename(file)), self.workdir + fastqc_foldername + '/' +
                        os.path.basename(file))

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nFASTQ UTILS' + '\033[0m')

        # Arguments
        self.args = args
        self.fastq_filenames = args.fastq_filenames
        self.fastqc = args.fastqc
        self.fasta = args.fasta
        self.qualitytrim = args.qualitytrim
        self.dereplicate = args.dereplicate
        self.merge = args.merge

        # Metadata
        self.num_reads = len(self.fastq_filenames)
        self.workdir = os.path.dirname(self.fastq_filenames[0])
        self.metagenome_id = self.fastq_filenames[0][:11]

        # BBDuk directory grab
        cmd = 'which bbduk.sh'
        self.bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
        self.bbduk_dir = self.bbduk_dir.split('/')[:-1]
        self.bbduk_dir = '/'.join(self.bbduk_dir)

        # User input validation
        if self.num_reads > 2:
            print('\nInvalid number of arguments for .fastq files passed. '
                  'Exiting.')
            quit()

        if self.dereplicate and self.qualitytrim is False:
            print('\nCannot run dereplication unless quality trimming with the -qt flag is specified as well. '
                  'Exiting.')
            quit()

        # Run methods
        # -qt -dr -fc
        if self.qualitytrim and self.fastqc and self.dereplicate:
            print('\033[92m' + '\033[1m' + '\nRunning BBDuk ==> Dereplication ==> FastQC pipeline... ' + '\033[0m')

            if self.num_reads == 2:
                read1_filtered, read2_filtered = self.quality_trim_bbduk()
                read1_deduped, read2_deduped = self.run_dedupe(read1_filtered, read2_filtered)
                self.run_fastqc(read1=read1_deduped, read2=read2_deduped)
            elif self.num_reads == 1:
                read1_filtered = self.quality_trim_bbduk()
                read1_deduped = self.run_dedupe(read1=read1_filtered)
                self.run_fastqc(read1=read1_deduped)

        # -qt -fc
        elif self.qualitytrim and self.fastqc:
            print('\033[92m' + '\033[1m' + '\nRunning BBDuk ==> FastQC pipeline...' + '\033[0m')
            read1_filtered, read2_filtered = self.quality_trim_bbduk()
            self.run_fastqc(read1=read1_filtered, read2=read2_filtered)

        # -qt
        elif self.qualitytrim:
            self.quality_trim_bbduk()

        # -fc
        elif self.fastqc:
            self.run_fastqc()

        # -fa
        if self.fasta:
            self.fasta()

        if self.merge and self.num_reads == 2:
            run_merge(self.fastq_filenames[0], self.fastq_filenames[1])

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_filenames',
                        help='Path to the FastQ file(s) you would like to perform QC functions on. '
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
    parser.add_argument('-m', '--merge',
                        default=False,
                        action='store_true',
                        help='Runs bbmerge.sh against a pair of reads. '
                             'Will merge a read pair into a single .fastq.gz file. '
                             'Run this without any other flags.')
    arguments = parser.parse_args()

    x = FastQUtils(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    print('\033[92m' + '\033[1m' + '\nFinished FastQUtils functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')

