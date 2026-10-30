#!/usr/bin/env python3

"""
In silico Bruce-ladder, Suis-ladder and AMOS PCR, run end to end from assembled genomes.

For every FASTA file the primer pairs of each requested assay are located in the sequence, the
resulting amplicon sizes are read as a gel would be read, and a genotype is called.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from insilico_pcr import InSilicoPCR, Amplicon
from genotype import Genotyper
import multiprocessing
import assays
import sys
import os

__author__ = 'duceppemo'
__version__ = '1.0'


def process_sample(job):
    """
    Run every requested assay against one FASTA file. Module level so it can be sent to a worker.
    :param job: type TUPLE: (FASTA file, list of assay names, InSilicoPCR object)
    :return: type TUPLE: (sample name, list of Amplicon named tuples)
    """
    fasta_file, assay_list, pcr = job
    amplicons = list()
    for assay in assay_list:
        amplicons.extend(pcr.amplify(fasta_file=fasta_file, assay=assay, primers=assays.PRIMERS[assay]))
    return InSilicoPCR.strain_namer(fasta_file=fasta_file), amplicons


class BruceLadder(object):

    def __init__(self, args):
        # Command line arguments
        self.input = args.input
        self.outdir = args.outdir
        self.threads = args.threads
        self.tolerance = args.tolerance

        # Default to every assay when none was requested
        self.assays = [assay for assay, requested in (('bruce', args.bruce),
                                                      ('suis', args.suis),
                                                      ('amos', args.amos)) if requested]
        if not self.assays:
            self.assays = assays.assay_names()

        self.pcr = InSilicoPCR(mismatches=args.mismatches,
                               three_prime=args.three_prime,
                               min_size=args.min_size,
                               max_size=args.max_size)

        # Data objects
        self.fasta_files = list()
        self.amplicons = list()
        self.calls = dict()

        # Run script
        self.run()

    def run(self):
        self.find_input()
        self.run_pcr()
        self.call_genotypes()
        self.write_reports()

    def find_input(self):
        """
        Locate the FASTA files to process
        """
        if not os.path.exists(self.input):
            raise FileNotFoundError('Cannot find the input path: {path}'.format(path=self.input))
        self.fasta_files = InSilicoPCR.file_list(path=self.input)
        if not self.fasta_files:
            raise FileNotFoundError('No FASTA files found in: {path}'.format(path=self.input))
        print('Found {n} FASTA file(s)'.format(n=len(self.fasta_files)), file=sys.stderr)

    def run_pcr(self):
        """
        Perform the in silico PCR, one worker per FASTA file
        """
        print('Running in silico PCR for assay(s): {assays}'.format(assays=', '.join(self.assays)),
              file=sys.stderr)
        jobs = [(fasta_file, self.assays, self.pcr) for fasta_file in self.fasta_files]
        threads = min(self.threads, len(jobs))
        if threads > 1:
            with multiprocessing.Pool(processes=threads) as pool:
                results = pool.map(process_sample, jobs)
        else:
            results = [process_sample(job) for job in jobs]
        for sample, amplicons in results:
            self.amplicons.extend(amplicons)
            print('  {sample}: {n} amplicon(s)'.format(sample=sample, n=len(amplicons)),
                  file=sys.stderr)

    def call_genotypes(self):
        """
        Convert the amplicon sizes of each sample into a genotype call, per assay
        """
        print('Calling genotypes', file=sys.stderr)
        samples = [InSilicoPCR.strain_namer(fasta_file=f) for f in self.fasta_files]
        for assay in self.assays:
            genotyper = Genotyper(assay=assay, tolerance=self.tolerance)
            for sample in samples:
                # Distinct sizes only: a locus repeated in the genome yields several products of
                # the same size, which a gel shows as a single band
                sizes = sorted({amplicon.size for amplicon in self.amplicons
                                if amplicon.sample == sample and amplicon.assay == assay},
                               reverse=True)
                call, pattern, extra = genotyper.call(sizes=sizes)
                self.calls[(sample, assay)] = {'call': call,
                                               'sizes': sizes,
                                               'pattern': pattern,
                                               'extra': extra,
                                               'bands': genotyper.bands}

    def write_reports(self):
        """
        Write the amplicon table, the consolidated call table, and one two column report per assay.
        The two column reports keep the layout the bruce_rates.py and suis_rates.py scripts expect.
        """
        os.makedirs(self.outdir, exist_ok=True)
        amplicon_report = os.path.join(self.outdir, 'amplicons.tsv')
        with open(amplicon_report, 'w') as out:
            out.write('\t'.join(Amplicon._fields) + '\n')
            for amplicon in sorted(self.amplicons, key=lambda a: (a.sample, a.assay, a.size)):
                out.write('\t'.join(str(field) for field in amplicon) + '\n')

        call_report = os.path.join(self.outdir, 'calls.tsv')
        with open(call_report, 'w') as out:
            out.write('\t'.join(['sample', 'assay', 'call', 'bands', 'sizes', 'extra_bands']) + '\n')
            for (sample, assay), result in sorted(self.calls.items()):
                bands = ';'.join('{low}-{high}'.format(low=low, high=high) if high > low
                                 else str(low)
                                 for (low, high), present in zip(result['bands'], result['pattern'])
                                 if present)
                out.write('\t'.join([sample,
                                     assay,
                                     result['call'],
                                     bands,
                                     ';'.join(str(size) for size in result['sizes']),
                                     ';'.join(str(size) for size in result['extra'])]) + '\n')

        for assay in self.assays:
            assay_report = os.path.join(self.outdir, 'report_{assay}.tsv'.format(assay=assay))
            with open(assay_report, 'w') as out:
                for (sample, sample_assay), result in sorted(self.calls.items()):
                    if sample_assay == assay:
                        out.write('{sample}\t{call}\n'.format(sample=sample, call=result['call']))
        print('Reports written to {outdir}'.format(outdir=os.path.abspath(self.outdir)), file=sys.stderr)


if __name__ == '__main__':
    parser = ArgumentParser(description='Perform in silico Bruce-ladder, Suis-ladder and/or AMOS PCR '
                                        'on assembled Brucella genomes',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='genomes/',
                        required=True,
                        help='Input FASTA file, or folder of FASTA files, gzipped or not. '
                             'One sample per file.')
    parser.add_argument('-o', '--outdir', metavar='results/',
                        required=True,
                        help='Output folder for the reports')
    parser.add_argument('--bruce', action='store_true',
                        help='Perform Bruce-ladder')
    parser.add_argument('--suis', action='store_true',
                        help='Perform Suis-ladder')
    parser.add_argument('--amos', action='store_true',
                        help='Perform AMOS PCR')
    parser.add_argument('-t', '--threads', metavar='N',
                        type=int, default=multiprocessing.cpu_count(),
                        help='Number of FASTA files to process in parallel')
    parser.add_argument('--mismatches', metavar='N',
                        type=int, default=2,
                        help='Mismatches tolerated per primer')
    parser.add_argument('--three-prime', metavar='N',
                        type=int, default=3,
                        help="Bases at the 3' end of a primer that must match exactly")
    parser.add_argument('--tolerance', metavar='BP',
                        type=int, default=5,
                        help='Size tolerance in bp when assigning an amplicon to an expected band')
    parser.add_argument('--min-size', metavar='BP',
                        type=int, default=50,
                        help='Smallest product to report, primers included')
    parser.add_argument('--max-size', metavar='BP',
                        type=int, default=3000,
                        help='Largest product to report, primers included')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    # Get the arguments into an object
    arguments = parser.parse_args()

    BruceLadder(arguments)
