#!/usr/bin/env python3

"""
Tests for the two legacy scoring scripts, bruce_rates.py and suis_rates.py.

Neither is a module: each reads a two column report, loops over a fixed species or biovar list and
prints a block of counts per name. These tests run them the way the Development page documents,
with the report as the first argument, and read the printed blocks back.

They exist because both scripts carried bugs that only the printed numbers reveal: every ratio was
divided without checking the denominator, and the counters that compare the call to the species
list did so with == against a list that holds substrings of the calls.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

__author__ = 'duceppemo'

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 'Total ID (pre-ID or by PCR) <name>: <n>' opens each block
HEADER = re.compile(r'^Total ID \(pre-ID or by PCR\) (?P<name>\S+): (?P<total>\d+)$')
# 'Genomes pre-ID as abort: 13', or 'Genomes pre-ID as suis bv01: 4' in suis_rates.py
PRE_ID = re.compile(r'^Genomes pre-ID as (?:suis )?\S+: (?P<pre_id>\d+)$')
# 'Wrong pre-ID: 2/11 (18.2%)', or '0/0' when there is nothing to divide by
RATIO = re.compile(r'^(?P<label>[^:]+): (?P<num>\d+)/(?P<den>\d+)(?: \((?P<pct>[\d.]+)%\))?$')


def run_rates(script, *args, **kwargs):
    """
    Run one of the rates scripts
    :param script: type STR: File name of the script, at the repository root
    :param args: Command line arguments to pass it
    :param kwargs: cwd, to run from somewhere other than the repository
    :return: type subprocess.CompletedProcess
    """
    command = [sys.executable, os.path.join(REPO, script)] + list(args)
    return subprocess.run(command, capture_output=True, text=True, check=True,
                          cwd=kwargs.get('cwd', REPO))


def parse(stdout):
    """
    Read the printed blocks back into a dictionary
    :param stdout: type STR: Everything the script printed
    :return: type DICT: name: {'total': anyID, 'pre_id': preID, 'ratios': {label: (num, den)}}
    """
    blocks = dict()
    current = None
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        header = HEADER.match(line)
        if header:
            current = {'total': int(header.group('total')), 'pre_id': None, 'ratios': dict()}
            blocks[header.group('name')] = current
            continue
        pre_id = PRE_ID.match(line)
        if pre_id and current is not None:
            current['pre_id'] = int(pre_id.group('pre_id'))
            continue
        ratio = RATIO.match(line)
        if ratio and current is not None:
            current['ratios'][ratio.group('label')] = (int(ratio.group('num')),
                                                       int(ratio.group('den')))
    return blocks


class RatesTestCase(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def report(self, rows, name='report.tsv'):
        """
        Write a two column report
        :param rows: type LIST: (file name, call) tuples, or raw strings to write verbatim
        :param name: type STR: File name to write
        :return: type STR: Path of the report
        """
        path = os.path.join(self.tmpdir, name)
        with open(path, 'w') as out:
            for row in rows:
                out.write(row if isinstance(row, str) else '{0}\t{1}\n'.format(*row))
        return path


class TestBruceRates(RatesTestCase):

    def test_a_report_of_agreeing_calls_scores_every_species(self):
        path = self.report([('abortus__A__1', 'abortus'),
                            ('melitensis__B__2', 'melitensis'),
                            ('ovis__C__3', 'ovis')])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        # One block per name in the script's species list
        self.assertEqual(sorted(blocks),
                         sorted(['abort', 'suis', 'ovis', 'canis', 'ceti', 'pedialis',
                                 'melitensis', 'neotomae']))
        for species in ['abort', 'melitensis', 'ovis']:
            self.assertEqual(blocks[species]['ratios']['ID by PCR'], (1, 1))
            self.assertEqual(blocks[species]['ratios']['Wrong pre-ID'], (0, 1))

    def test_a_call_naming_another_species_is_a_wrong_pre_id(self):
        # The counter that could never fire: the pre-ID says melitensis, the PCR says suis
        path = self.report([('melitensis__GTC_3P1517__1', 'suis_bv2/suis_bv1_3_4')])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        self.assertEqual(blocks['melitensis']['ratios']['Wrong pre-ID'], (1, 1))
        # and it is not also counted as unidentified
        self.assertEqual(blocks['melitensis']['ratios'][
            'False negatives (Pre-ID, but not ID by PCR, and not ID as other species)'], (0, 1))

    def test_a_substring_species_name_is_still_matched(self):
        # 'abort' has to be found inside 'abortus' and 'suis' inside 'suis_bv2/suis_bv1_3_4'.
        # Comparing the call to the species list with == misses both.
        path = self.report([('abortus__A__1', 'suis_bv2/suis_bv1_3_4')])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        self.assertEqual(blocks['abort']['ratios']['Wrong pre-ID'], (1, 1))
        self.assertEqual(blocks['suis']['ratios']['Just ID by PCR'], (1, 1))

    def test_an_unnamed_call_is_a_false_negative_not_a_wrong_pre_id(self):
        # N/A, RB51 and S19 name no field species, so they land in the other bucket
        for call in ['N/A', 'RB51', 'S19']:
            path = self.report([('abortus__A__1', call)], name='{0}.tsv'.format(call.replace('/', '')))
            blocks = parse(run_rates('bruce_rates.py', path).stdout)
            self.assertEqual(blocks['abort']['ratios']['Wrong pre-ID'], (0, 1),
                             'call {0}'.format(call))
            self.assertEqual(blocks['abort']['ratios'][
                'False negatives (Pre-ID, but not ID by PCR, and not ID as other species)'],
                (1, 1), 'call {0}'.format(call))

    def test_the_three_buckets_are_exclusive_and_exhaustive(self):
        # Every pre-identified genome lands in exactly one of three buckets: the PCR agrees, it
        # names another species, or it names none. So preID == agreed + wrong pre-ID + not
        # identified.
        #
        # Only checkable from stdout where anyID equals preID, because 'ID by PCR' and 'False
        # negatives' are both printed over anyID while 'Wrong pre-ID' is printed over preID. A
        # species named by a call on some other genome's row has anyID > preID and is skipped.
        path = self.report([('abortus__A__1', 'abortus'),
                            ('abortus__B__2', 'melitensis'),
                            ('abortus__C__3', 'N/A'),
                            ('ovis__D__4', 'ovis'),
                            ('canis__E__5', 'N/A')])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        checked = list()
        for species, block in blocks.items():
            ratios = block['ratios']
            wrong = ratios['Wrong pre-ID'][0]
            missed = ratios[
                'False negatives (Pre-ID, but not ID by PCR, and not ID as other species)'][0]
            agreed = ratios['ID by PCR'][0]
            if block['total'] != block['pre_id']:
                continue
            self.assertEqual(block['pre_id'], agreed + wrong + missed,
                             'species {0}: {1}'.format(species, ratios))
            checked.append((species, wrong))
        # Not vacuous: abort was checked and its wrong pre-ID bucket is the non-empty one
        self.assertIn(('abort', 1), checked)
        self.assertGreaterEqual(len(checked), 3)

    def test_overlapping_species_names_do_not_double_count(self):
        # 'ceti' and 'pedialis' both occur in the pinnipedialis_ceti call, so a correct call must
        # not also register as a wrong pre-ID for the other name
        path = self.report([('ceti__A__1', 'pinnipedialis_ceti'),
                            ('pinnipedialis__B__2', 'pinnipedialis_ceti')])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        self.assertEqual(blocks['ceti']['ratios']['Wrong pre-ID'], (0, 1))
        self.assertEqual(blocks['pedialis']['ratios']['Wrong pre-ID'], (0, 1))

    def test_a_species_absent_from_the_report_does_not_divide_by_zero(self):
        # 'abort' is the first name in the list, so a report without an abortus genome used to
        # kill the script on its first iteration
        path = self.report([('ovis__A__1', 'ovis'), ('ovis__B__2', 'ovis')])
        result = run_rates('bruce_rates.py', path)
        blocks = parse(result.stdout)
        self.assertEqual(blocks['abort']['total'], 0)
        self.assertEqual(blocks['abort']['ratios']['ID by PCR'], (0, 0))
        # The percentage is dropped rather than divided
        self.assertIn('ID by PCR: 0/0', result.stdout)
        self.assertNotIn('nan', result.stdout.lower())

    def test_an_empty_report_scores_every_species_at_zero(self):
        path = self.report([])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        self.assertEqual(len(blocks), 8)
        for species, block in blocks.items():
            self.assertEqual(block['total'], 0, species)

    def test_blank_lines_are_skipped(self):
        path = self.report(['\n', 'abortus__A__1\tabortus\n', '\n', '   \n'])
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        self.assertEqual(blocks['abort']['total'], 1)

    def test_a_report_with_no_trailing_newline_is_read(self):
        path = os.path.join(self.tmpdir, 'no_newline.tsv')
        with open(path, 'w') as out:
            out.write('abortus__A__1\tabortus')
        blocks = parse(run_rates('bruce_rates.py', path).stdout)
        self.assertEqual(blocks['abort']['total'], 1)

    def test_the_default_report_name_is_used_when_no_argument_is_given(self):
        self.report([('abortus__A__1', 'abortus')], name='report_bruce.tsv')
        blocks = parse(run_rates('bruce_rates.py', cwd=self.tmpdir).stdout)
        self.assertEqual(blocks['abort']['total'], 1)


class TestSuisRates(RatesTestCase):

    def test_each_biovar_is_scored(self):
        path = self.report([('suis_bv-1__A__1', 'bv01'),
                            ('suis_bv-2__B__2', 'bv02'),
                            ('suis_bv-5__C__3', 'bv05')])
        blocks = parse(run_rates('suis_rates.py', path).stdout)
        self.assertEqual(sorted(blocks), ['bv01', 'bv02', 'bv03', 'bv04', 'bv05'])
        for biovar in ['bv01', 'bv02', 'bv05']:
            self.assertEqual(blocks[biovar]['ratios']['ID by PCR'], (1, 1))
            self.assertEqual(blocks[biovar]['ratios']['False negatives'], (0, 1))

    def test_a_disagreeing_biovar_is_a_false_negative(self):
        path = self.report([('suis_bv-1__A__1', 'bv03')])
        blocks = parse(run_rates('suis_rates.py', path).stdout)
        self.assertEqual(blocks['bv01']['ratios']['False negatives'], (1, 1))
        # and bv03 gained a call with no suis pre-ID behind it
        self.assertEqual(blocks['bv03']['ratios']['Just ID by PCR'], (0, 1))

    def test_a_call_without_a_suis_pre_id_is_a_false_positive(self):
        path = self.report([('suis_bv-2__A__1', 'bv02'), ('melitensis__B__2', 'bv02')])
        blocks = parse(run_rates('suis_rates.py', path).stdout)
        self.assertEqual(blocks['bv02']['ratios']['False positives'], (1, 1))
        self.assertEqual(blocks['bv02']['ratios']['Just ID by PCR'], (1, 2))

    def test_a_biovar_absent_from_the_report_does_not_divide_by_zero(self):
        path = self.report([('suis_bv-1__A__1', 'bv01')])
        result = run_rates('suis_rates.py', path)
        blocks = parse(result.stdout)
        self.assertEqual(blocks['bv04']['total'], 0)
        self.assertEqual(blocks['bv04']['ratios']['ID by PCR'], (0, 0))
        # With no pre-ID there is nothing to report a rate against, so those lines are omitted
        self.assertNotIn('False positives', blocks['bv04']['ratios'])

    def test_an_empty_report_scores_every_biovar_at_zero(self):
        path = self.report([])
        blocks = parse(run_rates('suis_rates.py', path).stdout)
        self.assertEqual(len(blocks), 5)
        for biovar, block in blocks.items():
            self.assertEqual(block['total'], 0, biovar)

    def test_blank_lines_are_skipped(self):
        path = self.report(['\n', 'suis_bv-1__A__1\tbv01\n', '\n'])
        blocks = parse(run_rates('suis_rates.py', path).stdout)
        self.assertEqual(blocks['bv01']['total'], 1)

    def test_the_default_report_name_is_used_when_no_argument_is_given(self):
        self.report([('suis_bv-1__A__1', 'bv01')], name='report_suis.tsv')
        blocks = parse(run_rates('suis_rates.py', cwd=self.tmpdir).stdout)
        self.assertEqual(blocks['bv01']['total'], 1)


if __name__ == '__main__':
    unittest.main()
