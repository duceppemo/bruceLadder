#!/usr/bin/env python3

"""
End to end tests: synthetic genomes in, reports out.

Every genotype of every assay is built as a synthetic genome named after the genotype, the CLI is
run over the folder, and the reported call is checked against the name. This exercises the whole
chain, it is not a validation of the assays against real Brucella sequence.
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from insilico_pcr import Amplicon
from tests.synthetic import build_genome, write_fasta, write_profile_genomes
import assays

__author__ = 'duceppemo'

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, 'bruce_ladder.py')

# Genotypes a pattern cannot separate, reported as the matching profiles joined by '/'
MERGED = {'suis_bv2': 'suis_bv1_3_4/suis_bv2',
          'suis_bv1_3_4': 'suis_bv1_3_4/suis_bv2'}


def run_cli(indir, outdir, *flags):
    """
    Run the CLI over a folder of FASTA files
    :param indir: type STR: Folder holding the FASTA files
    :param outdir: type STR: Folder for the reports
    :param flags: Extra command line arguments
    :return: type subprocess.CompletedProcess
    """
    command = [sys.executable, SCRIPT, '-i', indir, '-o', outdir] + list(flags)
    return subprocess.run(command, capture_output=True, text=True, check=True)


def read_calls(outdir):
    """
    Read calls.tsv into a dictionary
    :param outdir: type STR: Folder holding the reports
    :return: type DICT: (sample, assay): row dictionary
    """
    with open(os.path.join(outdir, 'calls.tsv')) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        calls = dict()
        for line in handle:
            row = dict(zip(header, line.rstrip('\n').split('\t')))
            calls[(row['sample'], row['assay'])] = row
        return calls


def normalised(call):
    """
    Sort the genotypes of a '/' joined call so that the test does not depend on profile order
    :param call: type STR: Call as reported
    :return: type STR: Genotypes sorted and rejoined
    """
    return '/'.join(sorted(call.split('/')))


class AssayCallsMixin(object):
    """
    Build one genome per genotype of self.assay and check every call
    """

    assay = None

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.indir = os.path.join(self.tmpdir, 'genomes')
        self.outdir = os.path.join(self.tmpdir, 'results')
        self.genotypes = write_profile_genomes(assay=self.assay, outdir=self.indir)
        run_cli(self.indir, self.outdir, '--{assay}'.format(assay=self.assay), '-t', '1')
        self.calls = read_calls(self.outdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_one_call_per_genome(self):
        self.assertEqual(sorted(sample for sample, _ in self.calls), self.genotypes)

    def test_every_genotype_is_called_correctly(self):
        for genotype in self.genotypes:
            expected = normalised(MERGED.get(genotype, genotype))
            self.assertEqual(normalised(self.calls[(genotype, self.assay)]['call']), expected,
                             'genome built as {genotype}'.format(genotype=genotype))

    def test_no_unexpected_amplicons(self):
        # A primer pair matching a locus it was not meant to amplify would show up here
        for genotype in self.genotypes:
            self.assertEqual(self.calls[(genotype, self.assay)]['extra_bands'], '',
                             'genome built as {genotype}'.format(genotype=genotype))

    def test_reported_sizes_are_the_profile_sizes(self):
        for genotype in self.genotypes:
            sizes = [int(size) for size
                     in self.calls[(genotype, self.assay)]['sizes'].split(';')]
            self.assertEqual(sorted(sizes), sorted(assays.PROFILES[self.assay][genotype]))


class TestBruceLadder(AssayCallsMixin, unittest.TestCase):
    assay = 'bruce'


class TestSuisLadder(AssayCallsMixin, unittest.TestCase):
    assay = 'suis'


class TestAmosPcr(AssayCallsMixin, unittest.TestCase):
    assay = 'amos'


class TestReports(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.indir = os.path.join(self.tmpdir, 'genomes')
        self.outdir = os.path.join(self.tmpdir, 'results')
        os.makedirs(self.indir)
        write_fasta(build_genome(assay='bruce', sizes=assays.PROFILES['bruce']['melitensis']),
                    os.path.join(self.indir, 'melitensis.fasta'))

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_amplicon_report_header_matches_the_named_tuple(self):
        run_cli(self.indir, self.outdir, '--bruce', '-t', '1')
        with open(os.path.join(self.outdir, 'amplicons.tsv')) as handle:
            self.assertEqual(handle.readline().rstrip('\n').split('\t'), list(Amplicon._fields))

    def test_amplicon_report_has_one_row_per_product(self):
        run_cli(self.indir, self.outdir, '--bruce', '-t', '1')
        with open(os.path.join(self.outdir, 'amplicons.tsv')) as handle:
            rows = handle.read().splitlines()[1:]
        self.assertEqual(len(rows), len(assays.PROFILES['bruce']['melitensis']))

    def test_per_assay_report_keeps_the_two_column_layout(self):
        # bruce_rates.py and suis_rates.py parse 'sample<TAB>call' with no header
        run_cli(self.indir, self.outdir, '--bruce', '-t', '1')
        with open(os.path.join(self.outdir, 'report_bruce.tsv')) as handle:
            rows = [line.rstrip('\n').split('\t') for line in handle if line.strip()]
        self.assertEqual(rows, [['melitensis', 'melitensis']])

    def test_default_runs_every_assay(self):
        run_cli(self.indir, self.outdir, '-t', '1')
        self.assertEqual(sorted(assay for _, assay in read_calls(self.outdir)),
                         sorted(assays.assay_names()))
        for assay in assays.assay_names():
            self.assertTrue(os.path.isfile(
                os.path.join(self.outdir, 'report_{assay}.tsv'.format(assay=assay))))

    def test_outdir_is_created(self):
        nested = os.path.join(self.outdir, 'deeper')
        run_cli(self.indir, nested, '--bruce', '-t', '1')
        self.assertTrue(os.path.isfile(os.path.join(nested, 'calls.tsv')))


class TestInputHandling(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.indir = os.path.join(self.tmpdir, 'genomes')
        self.outdir = os.path.join(self.tmpdir, 'results')
        os.makedirs(self.indir)
        self.contigs = build_genome(assay='suis', sizes=assays.PROFILES['suis']['bv03'])

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_gzipped_input(self):
        path = os.path.join(self.indir, 'bv03.fasta')
        write_fasta(self.contigs, path)
        with open(path, 'rb') as handle, gzip.open(path + '.gz', 'wb') as out:
            out.write(handle.read())
        os.remove(path)
        run_cli(self.indir, self.outdir, '--suis', '-t', '1')
        self.assertEqual(read_calls(self.outdir)[('bv03', 'suis')]['call'], 'bv03')

    def test_single_file_as_input(self):
        path = write_fasta(self.contigs, os.path.join(self.indir, 'bv03.fasta'))
        run_cli(path, self.outdir, '--suis', '-t', '1')
        self.assertEqual(read_calls(self.outdir)[('bv03', 'suis')]['call'], 'bv03')

    def test_several_workers_give_the_same_calls(self):
        for genotype, sizes in assays.PROFILES['suis'].items():
            write_fasta(build_genome(assay='suis', sizes=sizes),
                        os.path.join(self.indir, '{g}.fasta'.format(g=genotype)))
        serial = os.path.join(self.outdir, 'serial')
        parallel = os.path.join(self.outdir, 'parallel')
        run_cli(self.indir, serial, '--suis', '-t', '1')
        run_cli(self.indir, parallel, '--suis', '-t', '4')
        self.assertEqual(read_calls(serial), read_calls(parallel))

    def test_a_genome_without_any_product_is_not_a_call(self):
        write_fasta([('contig_1', 'ACGT' * 500)], os.path.join(self.indir, 'empty.fasta'))
        run_cli(self.indir, self.outdir, '--suis', '-t', '1')
        self.assertEqual(read_calls(self.outdir)[('empty', 'suis')]['call'], 'N/A')

    def test_missing_input_path_fails(self):
        with self.assertRaises(subprocess.CalledProcessError) as raised:
            run_cli(os.path.join(self.tmpdir, 'nope'), self.outdir, '--suis')
        self.assertIn('FileNotFoundError', raised.exception.stderr)

    def test_folder_without_fasta_files_fails(self):
        open(os.path.join(self.indir, 'notes.txt'), 'w').close()
        with self.assertRaises(subprocess.CalledProcessError) as raised:
            run_cli(self.indir, self.outdir, '--suis')
        self.assertIn('No FASTA files found', raised.exception.stderr)


if __name__ == '__main__':
    unittest.main()
