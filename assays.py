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
    # The three primers of Vemulapalli et al. 1999, Fig. 1, as two pairs sharing the forward
    # primer. Primer 3 is the shared forward primer, in wboA upstream of the insertion point.
    'rb51': {
        # Primer 3 + primer 1, both in wboA either side of the insertion point: 456 bp when wboA
        # is intact, 1300 bp in RB51, where the insertion adds 844 bp between them.
        'wboA': {'F': 'GCCAACCAACCCAAATGCTCACAA', 'R': 'TTAAGCGCTGATGCCATTTCCTTCAC'},
        # Primer 3 + primer 2. Primer 2 is chimeric: its 5' 22 bases are wboA and its 3' 12 bases
        # are IS711, so it anneals only where the two are fused. Each half occurs on its own in
        # the 75 genomes carrying the locus; the fused sequence exists in RB51 and nowhere else.
        'wboA_IS711': {'F': 'GCCAACCAACCCAAATGCTCACAA',
                       'R': 'TTTAGTTTGCCGTAATATAGGTCTAGAACCTGTC'},
    },
}


# Expected band sizes (bp) per genotype, carried over from the 2023 scripts except where a
# comment says otherwise. A genotype is called when every canonical band of the assay agrees
# in presence/absence with one of these profiles, within --tolerance bp.
#
# Order is load-bearing for any assay matched as a subset: Genotyper stops at the first profile
# that fits, so for amos and rb51 the listing order decides the call itself, not just the name.
# Profiles are listed most specific first so that the RB51/2308 profile is named ahead of the
# abortus profile it also satisfies, and S19, whose single band is a subset of both, is named
# last. For the exact-matched assays the order only affects reporting.
PROFILES = {
    'bruce': {
        'abortus': [1682, 774, 590, 450, 152],
        # The two abortus vaccine strains, each separated from abortus by a single band. RB51 is a
        # rough derivative of 2308 carrying IS711 inside wboA, which lengthens the 1682 bp product
        # of BMEI0998_BMEI0997 by the element: 2526 bp in the complete genome (GCF_011801185.1),
        # reported as ~2524 bp. S19 is missing the 590 bp band of BMEII0428, and keeps everything
        # else abortus has (GCF_000018725.1, GCF_038420195.1).
        'RB51': [2526, 774, 590, 450, 152],
        'S19': [1682, 774, 450, 152],
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
        # The RB51/2308 primer of Bricker and Halling (1995) lies inside IS711 itself, 711 bp
        # into the element, and pairs with the universal IS711 primer across the junction of two
        # tandem IS711 elements. The band therefore reports that arrangement and not RB51: strain
        # 2308 and its rough derivative RB51 give the same profile, and this assay cannot tell
        # them apart. The junction measures 387 bp in the complete genomes of both 2308
        # (GCF_000054005.1) and RB51 (GCF_011801185.1); the paper reports the product as 364 bp.
        # The 494 band belongs in this profile: B. ovis carries a tandem IS711 pair of its own
        # and yields the same 387 bp product, so the abortus band is what keeps it out.
        'RB51_or_2308': [494, 387, 178],
        'abortus_bv1_2_4': [494, 178],
        'melitensis': [733, 178],
        'ovis': [976, 178],
        'suis_bv1': [285, 178],
        # S19 carries a deletion that removes the eri target of the 178 bp band, so it is the one
        # genotype of this assay without it. Listed last because its single band is a subset of
        # the profiles above.
        'S19': [494],
    },
    # A confirmatory assay for one question only: is this RB51 or not. The 919 bp junction band is
    # diagnostic on its own and 1300 is not a canonical band of the assay, so it is tolerated as an
    # extra band under either match mode. What subset mode actually buys is the sample showing both
    # 456 and 919 - a mixed culture, or a co-assembly of RB51 with a smooth strain - which exact
    # matching would drop to N/A. RB51 is listed first so that such a sample is called RB51 rather
    # than not_RB51: do not reorder these two. A genome with neither band has no readable wboA
    # locus and is called N/A.
    #
    # Two limits of matching on size alone, neither reachable at sane settings, both recorded so
    # that a later edit does not walk into them. Profiles name sizes and never loci, so nothing
    # requires the 919 bp band to have come from the wboA_IS711 pair: an insertion of about 463 bp
    # at the same site would make the plain wboA pair yield 919 bp and call RB51. No such Brucella
    # element is known, and amplicons.tsv always records which locus gave the band. Separately,
    # the two bands lie 463 bp apart, so a --tolerance of 463 or more merges them into one and
    # every smooth genome calls RB51; the default is 5.
    'rb51': {
        'RB51': [919],
        'not_RB51': [456],
    },
}


# How observed bands are compared to a profile.
#   exact  : every canonical band of the assay must agree in presence/absence
#   subset : the bands of the profile must be present, additional bands are tolerated and the first
#            matching profile wins. Needed for AMOS, where RB51 and 2308 carry the abortus band
#            as well, and for rb51, where a sample carrying both bands must call RB51.
MATCH_MODE = {
    'bruce': 'exact',
    'suis': 'exact',
    'amos': 'subset',
    'rb51': 'subset',
}


# The Bruce-ladder band that reports the O-antigen locus. wboA and wboB sit inside the product of
# this pair, so a strain whose O-antigen locus is interrupted or deleted is rough and either shows
# the band at another size or loses it: RB51 carries IS711 inside wboA and runs at 2526 bp, while
# B. ovis and the deleted field isolate LBAB009 give no product at all. Only this one band is read;
# see the Validation page of the wiki for what that does and does not catch.
ROUGH_ASSAY = 'bruce'
ROUGH_LOCUS = 'BMEI0998_BMEI0997'
ROUGH_INTACT_SIZE = 1682


def assay_names():
    """
    :return: type LIST: Names of the available assays, in report order
    """
    return ['bruce', 'suis', 'amos', 'rb51']
