#!/usr/local/env python3

__author__ = 'duceppemo'
__version__ = '0.1'

from argparse import ArgumentParser


class SuisLadder(object):
    suis_pcr = {
        'BMEI1426_BMEI1427': {'F': 'TCGTCGGTGGACTGGATGAC', 'R': 'ATGGTCCGCAAGGTGCTTTT'},
        'BR1080': {'F': 'CCCTTGGTTTGTAGCGGTTG', 'R': 'TCATCGTCCTCCGTCATCCT'},
        'BMEI1688_MEI1687': {'F': 'TCAACTGCGTGAACAATGCT', 'R': 'GCGGGCTCTATCTCAAGGTC'},
        'BMEI0205': {'F': 'CGTCAACTCGCTGGCCAAGAG', 'R': 'GCAGGAGAACCGCAACCTAA'}
    }

    bruce_pcr = {
        'BMEI0998_BMEI0997': {'F': 'ATCCTATTGCCCCGATAAGG', 'R': 'GCTTCGCATTTTCACTGTAGC'},
        'BMEI0535_BMEI0536': {'F': 'GCGCATTCTTCGGTTATGAA', 'R': 'CGCAGGCGAAAACAGCTATAA'},
        'BMEII0843_BMEII0844': {'F': 'TTTACACAGGCAATCCAGCA', 'R': 'GCGTCCAGTTGTTGTTGATG'},
        'BMEI1426_BMEI1427': {'F': 'TCGTCGGTGGACTGGATGAC', 'R': 'ATGGTCCGCAAGGTGCTTTT'},
        'BMEII0428_BMEII0428': {'F': 'GCCGCTATTATGTGGACTGG', 'R': 'AATGACTTCACGGTCGTTCG'},
        'BR0953': {'F': 'GGAACACTACGCCACCTTGT', 'R': 'GATGGAGCAAACGCTGAAG'},
        'BMEI0752': {'F': 'CAGGCAAACCCTCAGAAGC', 'R': 'GATGTGGTAACGCACACCAA'},
        'BMEII0987': {'F': 'CGCAGACAGTGACCATCAAA', 'R': 'GTATTCAGCCCCCGTTACCT'}
    }

    amos_pcr = {
        'abortus_IS711': {'F': 'GACGAACGGAATTTTTCCAATCCC', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'melitensis_IS711': {'F': 'AAATCGCGTCCTTGCTGGTCTGA', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'ovis_IS711': {'F': 'CGGGTTCTGGCACCATCGTCG', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'suis_IS711': {'F': 'GCGCGGTTTTCTGAAGGTTCAGG', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'RB51_2308_IS711': {'F': 'CCCCGGAAGATATGCTTCGATCC', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'eri': {'F': 'GCGCCGCGAAGAACTTATCAA', 'R': 'CGCCATGTTAGCGGCGGTGA'}
    }

    # Band sizes expected for Bruce-ladder PCR
    bruce_results = {
        'abortus': [1682, 774, 590, 450, 152],
        'melitensis': [1682, 1071, 774, 590, 450, 152],
        'ovis': [1071, 774, 590, 450, 152],
        # 'microti': [1682, 1071, 774, 590, 450, 272, 152],
        'suis_bv2': [1682, 1071, 774, 590, 450, 272, 152],
        'suis_bv1_3_4': [1682, 1071, 774, 590, 451, 272, 152],
        'canis': [1682, 1071, 590, 451, 272, 152],
        'neotomae': [1682, 1071, 774, 590, 450, 272],
        'pinnipedialis_ceti': [1682, 1071, 774, 590, 152]  # sometimes one extra band of ~1100 to 1700bp
    }

    bruce_sizes = [1682, 1071, 774, 590, 450, 272, 152]
    bruce_size_range = list()
    for size in bruce_sizes:
        size_range = range(size - 5, size + 5)
        bruce_size_range.append(size_range)

    bruce_binary = {
        'abortus': [1, 0, 1, 1, 1, 0, 1],
        'melitensis': [1, 1, 1, 1, 1, 0, 1],
        'ovis': [0, 1, 1, 1, 1, 0, 1],
        'suis': [1, 1, 1, 1, 1, 1, 1],
        'canis': [1, 1, 0,  1, 1, 1, 1],
        'neotomae': [1, 1, 1, 1, 1, 0],
        'pinnipedialis_ceti': [1, 1, 1, 1, 0, 0, 1]  # sometimes one extra band of ~1100 to 1700bp
    }

    # Band sizes expected for Suis-ladder PCR
    suis_results = {
        'bv01': [774, 425, 197],
        'bv02': [774, 551, 278],
        'bv03': [774, 299, 197],
        'bv04': [774, 614, 197],
        'bv05': [774, 614, 278, 197]
    }

    amos_results = {
        'abortus_bv1_2_4': [494, 178],  # 494, 493, 1012
        'melitensis': [733, 178],  # 733, 732, 730
        'ovis': [976, 178],
        'suis_bv1': [285, 178],  # 276?
        'RB51': [364, 178]
    }
    #     'S19': [504]
    # }

    def __init__(self, args):
        # Command line arguments
        self.report = args.report
        self.suis = args.suis
        self.bruce = args.bruce
        self.amos = args.amos

        # Data objects
        self.pcr_dict = dict([])

        # Run script
        self.run()

    def run(self):
        SuisLadder.parse_report(self.report, self.pcr_dict)
        if self.suis:
            self.do_suis_ladder(self.pcr_dict)
        if self.bruce:
            self.do_bruce_ladder(self.pcr_dict)
        if self.amos:
            self.do_amos_pcr(self.pcr_dict)
        # if not self.suis or not self.bruce:
        #     raise Exception('Please choose "--suis" and/or "--bruce" flags')

    @staticmethod
    def parse_report(r, d):
        """
        Parse consolidated report from in silico PCR tool into dictionary
        :param r: report
        :param d: dictionary
        :return:
        """
        with open(r, 'r') as f:
            f.readline()  # skip header
            for line in f:
                line = line.rstrip()
                if line == "":
                    continue  # skip empty lines
                fields = line.split('\t')
                sample = fields[0]
                size = [int(fields[3])]

                # Add to dictionary
                try:
                    d[sample].extend(size)
                except KeyError:
                    d[sample] = size

    def do_suis_ladder(self, d):
        """
        Perform in silico suis ladder
        :param d: Populated dictionary
        :return:
        """
        output_file = self.report.split('.')[:-1][0] + '_' + 'suis' + '.tsv'

        with open(output_file, 'w') as f:
            for sample, size_list in d.items():
                size_list = sorted(size_list, reverse=True)
                found = 0
                for ident, size_ref in self.suis_results.items():
                    if size_list == size_ref:
                        found = 1
                        f. write('{}\t{}\n'.format(sample, ident))
                        break
                if found == 0:
                    f.write('{}\t{}\n'.format(sample, 'N/A'))

    def do_bruce_ladder(self, d):
        """
        Perform in silico bruce ladder, allowing 10bp size variation for PCR products
        :param d: Populated dictionary
        :return:
        """

        output_file = self.report.split('.')[:-1][0] + '_' + 'bruce' + '.tsv'

        with open(output_file, 'w') as f:
            for sample, size_list in d.items():
                sizes_dectected = list()
                size_list = sorted(size_list, reverse=True)

                for size_range in self.bruce_size_range:
                    found = [size in size_range for size in size_list]
                    if any(found):
                        sizes_dectected.append(1)
                    else:
                        sizes_dectected.append(0)

                # for size in size_list:
                #     found = [size in x for x in self.bruce_size_range]
                #     if any(found):
                #         sizes_dectected.append(1)
                #     else:
                #         sizes_dectected.append(0)

                # get key by value
                ident = [x for x, y in self.bruce_binary.items() if y == sizes_dectected]
                if ident:
                    f. write('{}\t{}\n'.format(sample,  ident[0]))
                else:
                    f.write('{}\t{}\n'.format(sample, 'N/A'))

    # def do_bruce_ladder(self, d):
    #     """
    #     Perform in silico bruce ladder
    #     :param d: Populated dictionary
    #     :return:
    #     """
    #
    #     output_file = self.report.split('.')[:-1][0] + '_' + 'bruce' + '.tsv'
    #
    #     with open(output_file, 'w') as f:
    #         for sample, size_list in d.items():
    #             size_list = sorted(size_list, reverse=True)
    #             found = 0
    #             for ident, size_ref in self.bruce_results.items():
    #                 if size_list == size_ref:
    #                     found = 1
    #                     f. write('{}\t{}\n'.format(sample, ident))
    #                     break
    #             if found == 0:
    #                 f.write('{}\t{}\n'.format(sample, 'N/A'))

    def do_amos_pcr(self, d):
        """

        :param d:
        :return:
        """
        output_file = self.report.split('.')[:-1][0] + '_' + 'amos' + '.tsv'

        with open(output_file, 'w') as f:
            for sample, size_list in d.items():
                size_list = sorted(size_list, reverse=False)

                for ident, size_ref in self.amos_results.items():
                    found = all([s in size_list for s in size_ref])
                    if found:
                        f.write('{}\t{}\n'.format(sample, ident))
                        break
                    else:
                        f.write('{}\t{}\n'.format(sample, '-'))

    def write_output(self, t, r):
        """

        :param t: analysis type (bruce or suis)
        :param r: result
        :return:
        """
        output_file = self.report.split('.')[:-1] + '_' + type + '.tsv'

        with open(output_file, 'w') as f:
            f.write('test\n')


if __name__ == '__main__':
    parser = ArgumentParser(description='Perform in silico Suis or Bruce Ladder using the consolidated report of'
                                        'in silico PCR tool')
    parser.add_argument('-r', '--report', metavar='consolidated_report.txt',
                        required=True,
                        help='Input folder with fastq file(s),gzipped or not')
    parser.add_argument('--suis', action="store_true",
                        help='Perform Suis Ladder')
    parser.add_argument('--bruce', action="store_true",
                        help='Perform Bruce Ladder')
    parser.add_argument('--amos', action="store_true",
                        help='Perform AMOS PCR')

    # Get the arguments into an object
    arguments = parser.parse_args()

    SuisLadder(arguments)
