import os
import subprocess
import time
import argparse
import collections
from FastQUtils import run_merge

# TODO: Remove bad hardcoded paths to metaphlan2.py, export2graphlan.py, graphlan_annotate.py, and graphlan.py


def generate_abundance_table(initial_profile_file, taxonomic_level):
    """
    Generates a tab delimited file representing relative abundances for a specific taxonomic level derived from
    the initial metaphlan2 profile output. Specify taxonomic level of interest according to the ordered dict below.
    """

    # Ordered dictionary setup
    taxonomic_levels = collections.OrderedDict()
    taxonomic_levels['kingdom'] = 'k__'
    taxonomic_levels['class'] = 'c__'
    taxonomic_levels['order'] = 'o__'
    taxonomic_levels['family'] = 'f__'
    taxonomic_levels['genus'] = 'g__'
    taxonomic_levels['species'] = 's__'
    taxonomic_levels['strain'] = 't__'

    # Ordered dict manipulation to figure out what the next key is... excessively verbose in an attempt at clarity
    desired_tax_rank = taxonomic_levels[taxonomic_level]
    index_tax_rank = tuple(taxonomic_levels).index(taxonomic_level)
    tax_key_list = list(taxonomic_levels.keys())
    next_tax_level = tax_key_list[index_tax_rank + 1]
    next_tax_rank = taxonomic_levels[next_tax_level]

    output_file = initial_profile_file.replace('_profile', '_profile_{}'.format(taxonomic_level))

    p = subprocess.Popen("grep -E '({0})| "
                         "(^ID)' {2} "
                         "| grep -v '{1}' "
                         "| sed 's/^.*{0}//g' > "
                         "{3}".format(desired_tax_rank, next_tax_rank, initial_profile_file, output_file),
                         shell=True,
                         executable='/bin/bash')
    p.wait()
    return output_file


class MetaPhlAn2(object):

    def run_metaphlan(self, single_read):
        print('\nRunning MetaPhlAn2...')
        output_filename = self.workdir+'/'+os.path.basename(single_read)[:10]+'_profile.txt'

        # This is very bad, but it works for now
        p = subprocess.Popen('/home/dussaultf/MetaPhlAn2/bin/python '  # run through the metaphlan2 venv
                             '/home/dussaultf/PycharmProjects/metaphlan2/metaphlan2.py '
                             '{} '
                             '--input_type fastq '
                            # '-t rel_ab_w_read_stats '  # provides rel. abundance but breaks GraPhlAn
                             '--nproc 8 > '  # 8 cores
                             '{}'.format(single_read, output_filename),
                             shell=True)
        p.wait()

        return output_filename

    @staticmethod
    def generate_species_abundance_table(overall_abundance_profile):
        print("\nGenerating species abundance table...")

        output_file = overall_abundance_profile.replace('_profile', '_species_profile')

        # Pull only the species level abundances from *_profile.txt for downstream cladogram
        p = subprocess.Popen("grep -E '(s__)| "
                             "(^ID)' {} "
                             "| grep -v 't__' "
                             "| sed 's/^.*s__//g' > "
                             "{}".format(overall_abundance_profile, output_file),
                             shell=True,
                             executable='/bin/bash')
        p.wait()
        return output_file

    def create_cladogram(self, overall_abundance_profile):

        # Step 1: Create GraPhlAn input files
        print("\nCreating GraPhlAn input files...")

        annotated_file = overall_abundance_profile.replace('profile', 'profile.annot')

        p = subprocess.Popen('python2 '
                             '/home/dussaultf/PycharmProjects/metaphlan2/utils/export2graphlan/export2graphlan.py '
                             '--skip_rows 1,2 '
                             '-i {} '
                             '--tree merged_abundance.tree.txt '
                             '--annotation {} '
                             '--most_abundant 100 '
                             '--abundance_threshold 1 '
                             '--least_biomarkers 10 '
                             '--annotations 5,6 '
                             '--external_annotations 7 '
                             '--min_clade_size 1'.format(overall_abundance_profile, annotated_file),
                             shell=True,
                             cwd=self.workdir,
                             executable='/bin/bash')
        p.wait()

        # Step 2: Create cladogram pieces
        print("\nCreating cladogram input files...")

        p = subprocess.Popen('python2 '
                             '/home/dussaultf/PycharmProjects/graphlan/graphlan_annotate.py '
                             '--annot {} '
                             'merged_abundance.tree.txt '
                             'merged_abundance.xml'.format(annotated_file),
                             shell=True,
                             cwd=self.workdir,
                             executable='/bin/bash')
        p.wait()

        # Step 3: Visualize cladogram
        print("\nVisualizing cladogram...")
        p = subprocess.Popen('python2 '
                             '/home/dussaultf/PycharmProjects/graphlan/graphlan.py '
                             '--dpi 300 '
                             'merged_abundance.xml '
                             'merged_abundance.png '
                             '--external_legends',
                             shell=True,
                             cwd=self.workdir,
                             executable='/bin/bash')
        p.wait()

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nMetaPhlAn2 + GraPhlAn PIPELINE' + '\033[0m')

        # Arguments
        self.args = args
        self.fastq_filenames = args.fastq_filenames

        # Metadata
        self.num_reads = len(self.fastq_filenames)
        self.workdir = os.path.dirname(self.fastq_filenames[0])

        # PEP-8 warning suppression
        cladogram_profile = ''

        # Figure out reads/merging
        if self.num_reads == 2:
            self.read1 = self.fastq_filenames[0]
            self.read2 = self.fastq_filenames[1]
            output_file = run_merge(self.read1, self.read2)
            cladogram_profile = self.run_metaphlan(output_file)
        elif self.num_reads == 1:
            self.single_read = self.fastq_filenames[0]
            cladogram_profile = self.run_metaphlan(self.single_read)
        else:
            print('\nInvalid number of reads entered. Quitting.')
            quit()

        # Run GraPhlAn cladogram components
        self.generate_species_abundance_table(cladogram_profile)
        self.create_cladogram(cladogram_profile)


if __name__ == '__main__':
    start = time.time()

    parser = argparse.ArgumentParser(description='This is a wrapper for MetaPhlAn2 and GraPhlAn. This script will take '
                                                 'single or paired reads and ultimately generate a species abundance '
                                                 'profile as well as a cladogram figure. If a read pair is submitted, '
                                                 'the reads will first be merged using BBmerge.')
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
