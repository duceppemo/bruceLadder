#!/usr/bin/env python3

"""
Primer sets and expected band profiles for the Brucella multiplex PCR assays.

This module is data only. Edit the profiles here to change how genotypes are
called; no other module hard-codes band sizes.
"""

__author__ = 'duceppemo'


# Primer pairs, keyed by locus name. Sequences are 5'->3' as ordered.
PRIMERS = {
    'suis': {
        'BMEI1426_BMEI1427': {'F': 'TCGTCGGTGGACTGGATGAC', 'R': 'ATGGTCCGCAAGGTGCTTTT'},
        'BR1080': {'F': 'CCCTTGGTTTGTAGCGGTTG', 'R': 'TCATCGTCCTCCGTCATCCT'},
        'BMEI1688_MEI1687': {'F': 'TCAACTGCGTGAACAATGCT', 'R': 'GCGGGCTCTATCTCAAGGTC'},
        'BMEI0205': {'F': 'CGTCAACTCGCTGGCCAAGAG', 'R': 'GCAGGAGAACCGCAACCTAA'},
    },
    'bruce': {
        'BMEI0998_BMEI0997': {'F': 'ATCCTATTGCCCCGATAAGG', 'R': 'GCTTCGCATTTTCACTGTAGC'},
        'BMEI0535_BMEI0536': {'F': 'GCGCATTCTTCGGTTATGAA', 'R': 'CGCAGGCGAAAACAGCTATAA'},
        'BMEII0843_BMEII0844': {'F': 'TTTACACAGGCAATCCAGCA', 'R': 'GCGTCCAGTTGTTGTTGATG'},
        'BMEI1426_BMEI1427': {'F': 'TCGTCGGTGGACTGGATGAC', 'R': 'ATGGTCCGCAAGGTGCTTTT'},
        'BMEII0428_BMEII0428': {'F': 'GCCGCTATTATGTGGACTGG', 'R': 'AATGACTTCACGGTCGTTCG'},
        'BR0953': {'F': 'GGAACACTACGCCACCTTGT', 'R': 'GATGGAGCAAACGCTGAAG'},
        'BMEI0752': {'F': 'CAGGCAAACCCTCAGAAGC', 'R': 'GATGTGGTAACGCACACCAA'},
        'BMEII0987': {'F': 'CGCAGACAGTGACCATCAAA', 'R': 'GTATTCAGCCCCCGTTACCT'},
    },
    'amos': {
        'abortus_IS711': {'F': 'GACGAACGGAATTTTTCCAATCCC', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'melitensis_IS711': {'F': 'AAATCGCGTCCTTGCTGGTCTGA', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'ovis_IS711': {'F': 'CGGGTTCTGGCACCATCGTCG', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'suis_IS711': {'F': 'GCGCGGTTTTCTGAAGGTTCAGG', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'RB51_2308_IS711': {'F': 'CCCCGGAAGATATGCTTCGATCC', 'R': 'TGCCGATCACTTAAGGGCCTTCAT'},
        'eri': {'F': 'GCGCCGCGAAGAACTTATCAA', 'R': 'CGCCATGTTAGCGGCGGTGA'},
    },
}


# Expected band sizes (bp) per genotype, carried over unchanged from the 2023
# scripts. A genotype is called when every canonical band of the assay agrees
# in presence/absence with one of these profiles, within --tolerance bp.
#
# Order matters only for reporting: profiles are listed most specific first so
# that RB51 is named ahead of the abortus profile it also satisfies.
PROFILES = {
    'bruce': {
        'abortus': [1682, 774, 590, 450, 152],
        'melitensis': [1682, 1071, 774, 590, 450, 152],
        'ovis': [1071, 774, 590, 450, 152],
        'suis_bv2': [1682, 1071, 774, 590, 450, 272, 152],
        'suis_bv1_3_4': [1682, 1071, 774, 590, 451, 272, 152],
        'canis': [1682, 1071, 590, 451, 272, 152],
        'neotomae': [1682, 1071, 774, 590, 450, 272],
        # Sometimes one extra band of ~1100 to 1700 bp.
        'pinnipedialis_ceti': [1682, 1071, 774, 590, 152],
        # 'microti' has the same profile as suis_bv2 and cannot be resolved by
        # this assay: [1682, 1071, 774, 590, 450, 272, 152]
    },
    'suis': {
        'bv01': [774, 425, 197],
        'bv02': [774, 551, 278],
        'bv03': [774, 299, 197],
        'bv04': [774, 614, 197],
        'bv05': [774, 614, 278, 197],
    },
    'amos': {
        'RB51': [364, 178],
        'abortus_bv1_2_4': [494, 178],
        'melitensis': [733, 178],
        'ovis': [976, 178],
        'suis_bv1': [285, 178],
        # 'S19': [504]
    },
}


# How observed bands are compared to a profile.
#   exact  : every canonical band of the assay must agree in presence/absence
#   subset : the bands of the profile must be present, additional bands are tolerated and the first
#            matching profile wins. Needed for AMOS, where RB51 carries the abortus band as well.
MATCH_MODE = {
    'bruce': 'exact',
    'suis': 'exact',
    'amos': 'subset',
}


def assay_names():
    """
    :return: type LIST: Names of the available assays, in report order
    """
    return ['bruce', 'suis', 'amos']
