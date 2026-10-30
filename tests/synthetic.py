#!/usr/bin/env python3

"""
Build synthetic genomes whose in silico PCR products reproduce a chosen band profile.

Each expected band size is placed on its own contig as a forward primer, a random spacer and the
reverse complement of the reverse primer, so that the contig yields exactly one product of that
size. Loci are assigned to sizes in a fixed order; which locus carries which band is arbitrary and
has no biological meaning. These genomes exercise the pipeline, they do not validate the assays
against real Brucella sequence.
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from insilico_pcr import InSilicoPCR
import assays

__author__ = 'duceppemo'

FLANK = 300


def spacer(rng, length):
    """
    :param rng: type random.Random: Seeded generator, so that genomes are reproducible
    :param length: type INT: Number of bases to generate
    :return: type STR: Random nucleotide sequence
    """
    return ''.join(rng.choice('ACGT') for _ in range(length))


def product_contig(rng, forward, reverse, size, flank=FLANK):
    """
    Lay out a single product of the requested size on a contig
    :param rng: type random.Random: Seeded generator
    :param forward: type STR: Forward primer, 5'->3'
    :param reverse: type STR: Reverse primer, 5'->3'
    :param size: type INT: Product size in bp, primers included
    :param flank: type INT: Bases of spacer either side of the product
    :return: type STR: Contig sequence
    """
    right = InSilicoPCR.reverse_complement(reverse)
    gap = size - len(forward) - len(right)
    if gap < 0:
        raise ValueError('A {size} bp product does not fit the primers it was asked to use'.format(
            size=size))
    return spacer(rng, flank) + forward + spacer(rng, gap) + right + spacer(rng, flank)


def build_genome(assay, sizes, seed=None):
    """
    Build the contigs of a genome yielding one product per requested size
    :param assay: type STR: Assay name, a key of assays.PRIMERS
    :param sizes: type LIST: Expected band sizes in bp
    :param seed: type INT: Seed for the spacer sequences, derived from the sizes when not given
    :return: type LIST: (contig name, sequence) tuples
    """
    primers = assays.PRIMERS[assay]
    loci = sorted(primers)
    if len(sizes) > len(loci):
        raise ValueError('The {assay} assay has {n} primer pairs, {m} bands were asked for'.format(
            assay=assay, n=len(loci), m=len(sizes)))
    rng = random.Random(seed if seed is not None else sum(sizes))
    contigs = list()
    for locus, size in zip(loci, sorted(sizes, reverse=True)):
        contigs.append(('{locus}_{size}'.format(locus=locus, size=size),
                        product_contig(rng, primers[locus]['F'], primers[locus]['R'], size)))
    return contigs


def write_fasta(contigs, fasta_file):
    """
    :param contigs: type LIST: (contig name, sequence) tuples
    :param fasta_file: type STR: Path of the FASTA file to write
    :return: type STR: fasta_file
    """
    with open(fasta_file, 'w') as out:
        for name, seq in contigs:
            out.write('>{name}\n{seq}\n'.format(name=name, seq=seq))
    return fasta_file


def write_profile_genomes(assay, outdir):
    """
    Write one genome per genotype of an assay, named after the genotype
    :param assay: type STR: Assay name
    :param outdir: type STR: Folder to write the FASTA files into
    :return: type LIST: Sorted genotype names written
    """
    os.makedirs(outdir, exist_ok=True)
    for genotype, sizes in assays.PROFILES[assay].items():
        write_fasta(build_genome(assay=assay, sizes=sizes),
                    os.path.join(outdir, '{genotype}.fasta'.format(genotype=genotype)))
    return sorted(assays.PROFILES[assay])
