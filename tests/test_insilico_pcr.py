#!/usr/bin/env python3

"""
Unit tests for the in silico PCR
"""

import gzip
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from insilico_pcr import InSilicoPCR
from tests.synthetic import product_contig, spacer, write_fasta

__author__ = 'duceppemo'

FORWARD = 'TCGTCGGTGGACTGGATGAC'
REVERSE = 'ATGGTCCGCAAGGTGCTTTT'


class TestHelpers(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual(InSilicoPCR.reverse_complement('ACGT'), 'ACGT')
        self.assertEqual(InSilicoPCR.reverse_complement('AAAC'), 'GTTT')
        self.assertEqual(InSilicoPCR.reverse_complement(FORWARD),
                         'GTCATCCAGTCCACCGACGA')

    def test_reverse_complement_ambiguity_codes(self):
        # Ambiguity codes have to survive a round trip, assemblies do carry them
        self.assertEqual(InSilicoPCR.reverse_complement(
            InSilicoPCR.reverse_complement('ACGTRYKMSWBDHVN')), 'ACGTRYKMSWBDHVN')

    def test_seed_blocks_cover_the_primer(self):
        for mismatches in range(0, 4):
            blocks = InSilicoPCR.seed_blocks(FORWARD, mismatches)
            self.assertEqual(len(blocks), mismatches + 1)
            # Blocks must be contiguous and reassemble the primer, otherwise a site carrying
            # mismatches in the gaps would never be seeded
            self.assertEqual(''.join(block for _, block in blocks), FORWARD)
            for offset, block in blocks:
                self.assertEqual(FORWARD[offset:offset + len(block)], block)
                self.assertTrue(block)

    def test_strain_namer(self):
        self.assertEqual(InSilicoPCR.strain_namer('/data/BRU-001.fasta'), 'BRU-001')
        self.assertEqual(InSilicoPCR.strain_namer('/data/BRU-001.fna.gz'), 'BRU-001')
        self.assertEqual(InSilicoPCR.strain_namer('BRU.001.fa'), 'BRU.001')


class TestFileList(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def touch(self, name):
        path = os.path.join(self.tmpdir, name)
        open(path, 'w').close()
        return path

    def test_finds_fasta_extensions_and_skips_others(self):
        for name in ['a.fasta', 'b.fa', 'c.fna', 'd.fas', 'e.seq', 'f.ffn', 'g.fasta.gz']:
            self.touch(name)
        for name in ['notes.txt', 'reads.fastq', 'reads.fastq.gz', 'table.tsv']:
            self.touch(name)
        found = [os.path.basename(f) for f in InSilicoPCR.file_list(path=self.tmpdir)]
        self.assertEqual(found, ['a.fasta', 'b.fa', 'c.fna', 'd.fas', 'e.seq', 'f.ffn',
                                 'g.fasta.gz'])

    def test_single_file_is_returned_as_is(self):
        path = self.touch('a.fasta')
        self.assertEqual(InSilicoPCR.file_list(path=path), [os.path.abspath(path)])


class TestFindSites(unittest.TestCase):

    def setUp(self):
        self.pcr = InSilicoPCR(mismatches=2, three_prime=3, min_size=50, max_size=3000)
        self.rng = random.Random(7)

    def sites(self, seq, primer, anchor_three_prime=True):
        return self.pcr.find_sites(seq=seq, primer=primer, anchor_three_prime=anchor_three_prime)

    def test_exact_site_is_found(self):
        seq = spacer(self.rng, 200) + FORWARD + spacer(self.rng, 200)
        self.assertEqual(self.sites(seq, FORWARD), [(200, 0)])

    def test_site_with_tolerated_mismatches_is_found(self):
        # Two mismatches well away from the 3' end
        mutated = 'A' + FORWARD[1:5] + 'A' + FORWARD[6:]
        seq = spacer(self.rng, 200) + mutated + spacer(self.rng, 200)
        self.assertEqual(self.sites(seq, FORWARD), [(200, 2)])

    def test_site_with_too_many_mismatches_is_rejected(self):
        mutated = 'AAA' + FORWARD[3:]
        seq = spacer(self.rng, 200) + mutated + spacer(self.rng, 200)
        self.assertEqual(self.sites(seq, FORWARD), [])

    def test_three_prime_mismatch_is_rejected(self):
        # A single mismatch, but inside the 3' bases that must match exactly
        mutated = FORWARD[:-1] + ('A' if FORWARD[-1] != 'A' else 'C')
        seq = spacer(self.rng, 200) + mutated + spacer(self.rng, 200)
        self.assertEqual(self.sites(seq, FORWARD), [])
        # The same mismatch is tolerated once the 3' end is no longer anchored there
        self.assertEqual(self.sites(seq, FORWARD, anchor_three_prime=False), [(200, 1)])

    def test_an_anchor_longer_than_the_primer_matches_it_in_full(self):
        # The anchor is clamped to the primer length, so asking for more bases than the primer has
        # can only ever anchor all of them. Unclamped, the slice wrapped and anchored fewer: a
        # mismatch at position 4 of this 20 bp primer was rejected at --three-prime 20 but
        # accepted at 21 through 30, making a larger anchor looser instead of stricter.
        mutated = FORWARD[:3] + ('A' if FORWARD[3] != 'A' else 'C') + FORWARD[4:]
        seq = spacer(self.rng, 200) + mutated + spacer(self.rng, 200)
        for three_prime in [len(FORWARD), len(FORWARD) + 1, len(FORWARD) + 10,
                            len(FORWARD) * 2]:
            pcr = InSilicoPCR(mismatches=2, three_prime=three_prime)
            self.assertEqual(pcr.find_sites(seq=seq, primer=FORWARD, anchor_three_prime=True), [],
                             'anchor of {n} bases let a mismatch through'.format(n=three_prime))
            self.assertEqual(pcr.find_sites(seq=seq, primer=FORWARD, anchor_three_prime=False), [],
                             'anchor of {n} bases let a mismatch through'.format(n=three_prime))

    def test_repeated_locus_yields_every_site(self):
        seq = (spacer(self.rng, 100) + FORWARD + spacer(self.rng, 100) + FORWARD +
               spacer(self.rng, 100))
        self.assertEqual([start for start, _ in self.sites(seq, FORWARD)], [100, 220])


class TestAmplifyContig(unittest.TestCase):

    def setUp(self):
        self.pcr = InSilicoPCR(mismatches=2, three_prime=3, min_size=50, max_size=3000)
        self.rng = random.Random(11)

    def amplify(self, seq):
        return self.pcr.amplify_contig(seq=seq, forward=FORWARD, reverse=REVERSE)

    def test_product_on_plus_strand(self):
        seq = product_contig(self.rng, FORWARD, REVERSE, 590, flank=100)
        products = self.amplify(seq)
        self.assertEqual(len(products), 1)
        strand, start, end, f_mismatches, r_mismatches = products[0]
        self.assertEqual(strand, '+')
        self.assertEqual((start, end), (100, 690))
        self.assertEqual(end - start, 590)
        self.assertEqual((f_mismatches, r_mismatches), (0, 0))

    def test_product_on_minus_strand(self):
        # The same locus sitting the other way round in the assembly
        seq = InSilicoPCR.reverse_complement(
            product_contig(self.rng, FORWARD, REVERSE, 590, flank=100))
        products = self.amplify(seq)
        self.assertEqual(len(products), 1)
        strand, start, end, f_mismatches, r_mismatches = products[0]
        self.assertEqual(strand, '-')
        self.assertEqual(end - start, 590)
        self.assertEqual((f_mismatches, r_mismatches), (0, 0))

    def test_mismatches_are_reported_against_the_right_primer(self):
        forward = 'A' + FORWARD[1:]
        seq = product_contig(self.rng, forward, REVERSE, 400, flank=100)
        products = self.amplify(seq)
        self.assertEqual(len(products), 1)
        _, _, _, f_mismatches, r_mismatches = products[0]
        self.assertEqual((f_mismatches, r_mismatches), (1, 0))

    def test_product_above_max_size_is_not_reported(self):
        pcr = InSilicoPCR(mismatches=2, three_prime=3, min_size=50, max_size=500)
        seq = product_contig(self.rng, FORWARD, REVERSE, 590, flank=100)
        self.assertEqual(pcr.amplify_contig(seq=seq, forward=FORWARD, reverse=REVERSE), [])

    def test_product_below_min_size_is_not_reported(self):
        pcr = InSilicoPCR(mismatches=2, three_prime=3, min_size=600, max_size=3000)
        seq = product_contig(self.rng, FORWARD, REVERSE, 590, flank=100)
        self.assertEqual(pcr.amplify_contig(seq=seq, forward=FORWARD, reverse=REVERSE), [])

    def test_primers_without_a_site_yield_nothing(self):
        self.assertEqual(self.amplify(spacer(self.rng, 2000)), [])


class TestReadFasta(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.rng = random.Random(13)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_plain_and_gzipped_fasta_read_the_same(self):
        contigs = [('contig_1', product_contig(self.rng, FORWARD, REVERSE, 590, flank=50))]
        plain = write_fasta(contigs, os.path.join(self.tmpdir, 'sample.fasta'))
        gzipped = os.path.join(self.tmpdir, 'sample.fasta.gz')
        with open(plain, 'rb') as handle, gzip.open(gzipped, 'wb') as out:
            out.write(handle.read())
        self.assertEqual(InSilicoPCR.read_fasta(plain), InSilicoPCR.read_fasta(gzipped))
        self.assertEqual(InSilicoPCR.read_fasta(plain), contigs)

    def test_sequence_is_upper_cased(self):
        # Assemblies often carry soft-masked, lower case sequence
        path = write_fasta([('contig_1', 'acgtACGT')], os.path.join(self.tmpdir, 'lower.fasta'))
        self.assertEqual(InSilicoPCR.read_fasta(path), [('contig_1', 'ACGTACGT')])


if __name__ == '__main__':
    unittest.main()
