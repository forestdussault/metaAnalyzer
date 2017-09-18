import os
import subprocess
import time
import argparse
from FastQUtils import run_merge

class MetaPhlAn2(object):

    def run_metaphlan(self, single_read):
        output = self.workdir+'/'+os.path.basename(single_read)[:10]+'_profile.txt'

        # This is very bad, but it works for now
        p = subprocess.Popen('/home/dussaultf/MetaPhlAn2/bin/python ' # run through the metaphlan2 venv
                             '/home/dussaultf/PycharmProjects/metaphlan2/metaphlan2.py '
                             '{} '
                             '--input_type fastq '
                             '--nproc 8 > ' # 8 cores
                             '{}'.format(single_read, output),
                             shell=True)
        p.wait()

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nMETAPHLAN2' + '\033[0m')

        # Arguments
        self.args = args
        self.fastq_filenames = args.fastq_filenames

        # Metadata
        self.num_reads = len(self.fastq_filenames)
        self.workdir = os.path.dirname(self.fastq_filenames[0])

        # Figure out reads/merging
        if self.num_reads == 2:
            self.read1 = self.fastq_filenames[0]
            self.read2 = self.fastq_filenames[1]
            output_file = run_merge(self.read1, self.read2)
            self.run_metaphlan(output_file)
        elif self.num_reads == 1:
            self.single_read = self.fastq_filenames[0]
            self.run_metaphlan(self.single_read)
        else:
            print('Invalid number of reads entered. Quitting.')
            quit()


if __name__ == '__main__':
    start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_filenames',
                        help='Path to the FastQ file(s) you would like to perform MetaPhlAn2 analysis on. '
                             'Enter a single read to begin the MetaPhlAn2 profiling process. '
                             'If two reads are passed they will be automatically merged beforehand.',
                        nargs='*')
    arguments = parser.parse_args()

    x = MetaPhlAn2(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    # Bold green time courtesy of Adam and Andrew
    print('\033[92m' + '\033[1m' + '\nFinished MetaPhlAn2 functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')