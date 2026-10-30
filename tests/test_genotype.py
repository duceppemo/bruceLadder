#!/usr/bin/env python3

"""
Unit tests for the genotype calling
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genotype import Genotyper
import assays

__author__ = 'duceppemo'


class TestCanonicalBands(unittest.TestCase):

    def test_sizes_within_tolerance_are_merged(self):
        bands = Genotyper.canonical_bands(profiles={'a': [100, 101], 'b': [200]}, tolerance=5)
        self.assertEqual(bands, [(100, 101), (200, 200)])

    def test_sizes_beyond_tolerance_stay_apart(self):
        bands = Genotyper.canonical_bands(profiles={'a': [100, 106], 'b': [200]}, tolerance=5)
        self.assertEqual(bands, [(100, 100), (106, 106), (200, 200)])

    def test_bruce_merges_the_450_451_pair(self):
        # The Bruce-ladder profiles differ by 1 bp between suis bv2 and bv1/3/4, which no gel
        # resolves, so those sizes have to end up in a single band
        bands = Genotyper(assay='bruce', tolerance=5).bands
        self.assertEqual(bands, [(152, 152), (272, 272), (450, 451), (590, 590), (774, 774),
                                 (1071, 1071), (1682, 1682), (2526, 2526)])

    def test_suis_bands_are_all_distinct(self):
        bands = Genotyper(assay='suis', tolerance=5).bands
        self.assertEqual(bands, [(197, 197), (278, 278), (299, 299), (425, 425), (551, 551),
                                 (614, 614), (774, 774)])


class TestBandAssignment(unittest.TestCase):

    def setUp(self):
        self.genotyper = Genotyper(assay='suis', tolerance=5)

    def test_size_on_the_nose(self):
        self.assertEqual(self.genotyper.bands[self.genotyper.band_index(size=774)], (774, 774))

    def test_size_at_the_edge_of_the_tolerance(self):
        self.assertIsNotNone(self.genotyper.band_index(size=769))
        self.assertIsNotNone(self.genotyper.band_index(size=779))

    def test_size_beyond_the_tolerance(self):
        self.assertIsNone(self.genotyper.band_index(size=768))
        self.assertIsNone(self.genotyper.band_index(size=780))

    def test_extra_bands_are_the_unassigned_sizes(self):
        self.assertEqual(self.genotyper.extra_bands(sizes=[774, 1500, 197, 900]), [900, 1500])

    def test_pattern_ignores_unassigned_sizes(self):
        self.assertEqual(self.genotyper.pattern(sizes=[774, 1500]),
                         self.genotyper.pattern(sizes=[774]))


class TestExactMatching(unittest.TestCase):
    """
    Bruce-ladder and Suis-ladder are called on the full presence/absence pattern
    """

    def test_every_suis_profile_calls_itself(self):
        genotyper = Genotyper(assay='suis', tolerance=5)
        for genotype, sizes in assays.PROFILES['suis'].items():
            call, _, extra = genotyper.call(sizes=sizes)
            self.assertEqual(call, genotype)
            self.assertEqual(extra, [])

    def test_resolvable_bruce_profiles_call_themselves(self):
        genotyper = Genotyper(assay='bruce', tolerance=5)
        for genotype in ['abortus', 'RB51', 'S19', 'melitensis', 'ovis', 'canis', 'neotomae',
                         'pinnipedialis_ceti']:
            call, _, _ = genotyper.call(sizes=assays.PROFILES['bruce'][genotype])
            self.assertEqual(call, genotype)

    def test_rb51_is_abortus_with_the_wboA_band_shifted(self):
        # IS711 sits inside wboA in RB51, lengthening the 1682 bp product to 2526. Everything else
        # is the abortus pattern, so that one band is what separates them.
        genotyper = Genotyper(assay='bruce', tolerance=5)
        abortus = assays.PROFILES['bruce']['abortus']
        rb51 = assays.PROFILES['bruce']['RB51']
        self.assertEqual(sorted(set(abortus) ^ set(rb51)), [1682, 2526])
        self.assertEqual(genotyper.call(sizes=rb51)[0], 'RB51')
        self.assertEqual(genotyper.call(sizes=abortus)[0], 'abortus')

    def test_s19_is_abortus_without_the_590_band(self):
        # S19 differs from abortus by the absence of the BMEII0428 band and nothing else
        genotyper = Genotyper(assay='bruce', tolerance=5)
        abortus = assays.PROFILES['bruce']['abortus']
        s19 = assays.PROFILES['bruce']['S19']
        self.assertEqual(sorted(set(abortus) - set(s19)), [590])
        self.assertEqual(sorted(set(s19) - set(abortus)), [])
        self.assertEqual(genotyper.call(sizes=s19)[0], 'S19')

    def test_the_vaccine_strains_do_not_shadow_abortus(self):
        # Each of the three patterns has to be distinct, or an abortus field strain could be
        # reported as a vaccine strain
        genotyper = Genotyper(assay='bruce', tolerance=5)
        patterns = {name: genotyper.expected[name] for name in ['abortus', 'RB51', 'S19']}
        self.assertEqual(len(set(patterns.values())), 3)

    def test_suis_biovars_are_reported_as_indistinguishable(self):
        # suis bv2 and bv1/3/4 differ only by the 450 vs 451 bp band, so Bruce-ladder can never
        # name one of them alone. Both profiles have to come back as the pair.
        genotyper = Genotyper(assay='bruce', tolerance=5)
        for genotype in ['suis_bv2', 'suis_bv1_3_4']:
            call, _, _ = genotyper.call(sizes=assays.PROFILES['bruce'][genotype])
            self.assertEqual(sorted(call.split('/')), ['suis_bv1_3_4', 'suis_bv2'])

    def test_an_extra_canonical_band_breaks_the_call(self):
        # abortus plus the 272 bp band matches no profile and must not be called abortus
        genotyper = Genotyper(assay='bruce', tolerance=5)
        call, _, extra = genotyper.call(sizes=assays.PROFILES['bruce']['abortus'] + [272])
        self.assertEqual(call, 'N/A')
        self.assertEqual(extra, [])

    def test_a_missing_band_breaks_the_call(self):
        genotyper = Genotyper(assay='suis', tolerance=5)
        call, _, _ = genotyper.call(sizes=[774, 425])
        self.assertEqual(call, 'N/A')

    def test_no_amplicons_is_not_a_call(self):
        for assay in assays.assay_names():
            call, _, extra = Genotyper(assay=assay, tolerance=5).call(sizes=[])
            self.assertEqual(call, 'N/A')
            self.assertEqual(extra, [])

    def test_sizes_within_tolerance_still_call(self):
        genotyper = Genotyper(assay='suis', tolerance=5)
        call, _, _ = genotyper.call(sizes=[772, 427, 195])
        self.assertEqual(call, 'bv01')

    def test_non_canonical_sizes_are_reported_but_do_not_break_the_call(self):
        genotyper = Genotyper(assay='suis', tolerance=5)
        call, _, extra = genotyper.call(sizes=[774, 425, 197, 1200])
        self.assertEqual(call, 'bv01')
        self.assertEqual(extra, [1200])


class TestSubsetMatching(unittest.TestCase):
    """
    AMOS is called on a subset: extra bands are tolerated and the first matching profile wins
    """

    def setUp(self):
        self.genotyper = Genotyper(assay='amos', tolerance=5)

    def test_every_amos_profile_calls_itself(self):
        for genotype, sizes in assays.PROFILES['amos'].items():
            call, _, _ = self.genotyper.call(sizes=sizes)
            self.assertEqual(call, genotype)

    def test_rb51_2308_wins_over_the_abortus_band_it_also_carries(self):
        # RB51 is derived from abortus 2308 and shows the abortus band as well, so profile order
        # in assays.PROFILES is what keeps it from being called abortus
        call, _, _ = self.genotyper.call(sizes=[494, 387, 178])
        self.assertEqual(call, 'RB51_or_2308')

    def test_abortus_without_the_tandem_is711_band_is_abortus(self):
        call, _, _ = self.genotyper.call(sizes=[494, 178])
        self.assertEqual(call, 'abortus_bv1_2_4')

    def test_rb51_2308_is_listed_before_abortus(self):
        order = list(assays.PROFILES['amos'])
        self.assertLess(order.index('RB51_or_2308'), order.index('abortus_bv1_2_4'))

    def test_s19_is_listed_last(self):
        # The S19 profile is the abortus band on its own, a subset of every other profile that
        # carries it, so anything but last place would shadow them
        order = list(assays.PROFILES['amos'])
        self.assertEqual(order.index('S19'), len(order) - 1)

    def test_the_tandem_is711_band_without_the_abortus_band_is_not_rb51(self):
        # B. ovis carries a tandem IS711 pair of its own and gives the same 387 bp product
        call, _, _ = self.genotyper.call(sizes=[976, 387, 178])
        self.assertEqual(call, 'ovis')

    def test_the_shared_band_alone_is_not_a_call(self):
        call, _, _ = self.genotyper.call(sizes=[178])
        self.assertEqual(call, 'N/A')

    def test_s19_is_the_only_amos_profile_without_the_178_band(self):
        # The 178 bp eri product is carried by every Brucella genome tested, which is what makes
        # its absence the S19 call: S19 carries a deletion that removes the target.
        for genotype, sizes in assays.PROFILES['amos'].items():
            if genotype == 'S19':
                self.assertNotIn(178, sizes)
            else:
                self.assertIn(178, sizes)


class TestAssayData(unittest.TestCase):

    def test_every_assay_has_primers_profiles_and_a_match_mode(self):
        for assay in assays.assay_names():
            self.assertIn(assay, assays.PRIMERS)
            self.assertIn(assay, assays.PROFILES)
            self.assertIn(assays.MATCH_MODE[assay], ['exact', 'subset'])

    def test_primers_are_unambiguous_dna(self):
        for assay, primers in assays.PRIMERS.items():
            for locus, pair in primers.items():
                for end in ['F', 'R']:
                    self.assertTrue(set(pair[end]) <= set('ACGT'),
                                    '{assay} {locus} {end} is not plain ACGT'.format(
                                        assay=assay, locus=locus, end=end))

    def test_no_assay_asks_for_more_bands_than_it_has_primer_pairs(self):
        # A profile with more bands than loci could never be produced
        for assay, profiles in assays.PROFILES.items():
            for genotype, sizes in profiles.items():
                self.assertLessEqual(len(sizes), len(assays.PRIMERS[assay]),
                                     '{assay} {genotype}'.format(assay=assay, genotype=genotype))

    def test_profiles_are_distinguishable_or_knowingly_merged(self):
        # Any two profiles sharing a pattern can only ever be reported together, so the set of
        # collisions is pinned here rather than discovered on real data
        collisions = set()
        for assay in assays.assay_names():
            genotyper = Genotyper(assay=assay, tolerance=5)
            seen = dict()
            for genotype, pattern in genotyper.expected.items():
                if pattern in seen:
                    collisions.add((assay, tuple(sorted([seen[pattern], genotype]))))
                seen[pattern] = genotype
        self.assertEqual(collisions,
                         {('bruce', ('suis_bv1_3_4', 'suis_bv2'))})


if __name__ == '__main__':
    unittest.main()
