#!/usr/bin/env python3

"""
Translate observed amplicon sizes into a genotype call.

The canonical bands of an assay are the sizes used by its profiles, with sizes that fall within the
size tolerance of each other merged into a single band: the Bruce-ladder profiles differ by 1 bp
between suis biovar 2 and biovars 1, 3 and 4, which no gel and no size tolerance can separate.
"""

import assays

__author__ = 'duceppemo'


def rough_flag(assay, wboa_sizes, tolerance=5):
    """
    Read the Bruce-ladder O-antigen band as a smooth/rough flag
    :param assay: type STR: Assay the sizes came from; only assays.ROUGH_ASSAY carries the band
    :param wboa_sizes: type LIST: Sizes observed for assays.ROUGH_LOCUS in one sample
    :param tolerance: type INT: Size tolerance in bp
    :return: type STR: 'no' when the band is intact, 'yes' when it is shifted or missing, and an
    empty string for an assay that does not carry the band
    """
    if assay != assays.ROUGH_ASSAY:
        return ''
    intact = any(abs(size - assays.ROUGH_INTACT_SIZE) <= tolerance for size in wboa_sizes)
    return 'no' if intact else 'yes'


class Genotyper(object):

    def __init__(self, assay, tolerance=5):
        """
        :param assay: type STR: Assay name, a key of assays.PROFILES
        :param tolerance: type INT: Size tolerance in bp when assigning an amplicon to a band
        """
        self.assay = assay
        self.tolerance = tolerance
        self.profiles = assays.PROFILES[assay]
        self.mode = assays.MATCH_MODE[assay]
        self.bands = self.canonical_bands(profiles=self.profiles, tolerance=tolerance)
        # Presence/absence of every canonical band, per profile
        self.expected = {name: self.pattern(sizes=sizes) for name, sizes in self.profiles.items()}

    @staticmethod
    def canonical_bands(profiles, tolerance):
        """
        Merge the sizes used by the profiles of an assay into bands
        :param profiles: type DICT: Dictionary of genotype name: list of expected sizes
        :param tolerance: type INT: Size tolerance in bp
        :return: type LIST: (lowest size, highest size) of each band, smallest band first
        """
        sizes = sorted({size for profile_sizes in profiles.values() for size in profile_sizes})
        bands = list()
        for size in sizes:
            if bands and size - bands[-1][1] <= tolerance:
                bands[-1][1] = size
            else:
                bands.append([size, size])
        return [(low, high) for low, high in bands]

    def band_index(self, size):
        """
        Find the canonical band an amplicon size belongs to
        :param size: type INT: Amplicon size in bp
        :return: type INT: Index into self.bands, or None when the size matches no canonical band
        """
        for index, (low, high) in enumerate(self.bands):
            if low - self.tolerance <= size <= high + self.tolerance:
                return index
        return None

    def pattern(self, sizes):
        """
        Convert a list of sizes into the presence/absence of each canonical band
        :param sizes: type LIST: Amplicon sizes in bp
        :return: type TUPLE: True/False per canonical band
        """
        present = [False] * len(self.bands)
        for size in sizes:
            index = self.band_index(size=size)
            if index is not None:
                present[index] = True
        return tuple(present)

    def extra_bands(self, sizes):
        """
        Report the sizes that belong to no canonical band of the assay
        :param sizes: type LIST: Amplicon sizes in bp
        :return: type LIST: Sorted sizes that were not assigned to a band
        """
        return sorted(size for size in sizes if self.band_index(size=size) is None)

    def call(self, sizes):
        """
        Name the genotype matching a set of observed amplicon sizes
        :param sizes: type LIST: Amplicon sizes in bp
        :return: type TUPLE: (call, observed pattern, unassigned sizes). The call is 'N/A' when no
        profile matches, and the matching profiles joined by '/' when the assay cannot separate them
        """
        observed = self.pattern(sizes=sizes)
        matches = list()
        for name, expected in self.expected.items():
            if self.mode == 'subset':
                # Only the bands the profile calls for have to be present
                if all(not want or seen for want, seen in zip(expected, observed)):
                    matches.append(name)
                    break
            elif observed == expected:
                matches.append(name)
        call = '/'.join(matches) if matches else 'N/A'
        return call, observed, self.extra_bands(sizes=sizes)
