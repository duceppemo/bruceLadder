#!/usr/bin/env python3

"""
In silico PCR on assembled sequences.

Primer sites are located directly in the nucleotide sequence, allowing a
configurable number of mismatches while keeping the 3' end of each primer
exact. Any forward site paired with a downstream reverse site on the same
sequence, within the allowed product size range, is reported as an amplicon.
"""

from collections import namedtuple
from Bio import SeqIO
import gzip
import os

__author__ = 'duceppemo'

# Reported coordinates are 1-based and inclusive; size includes both primers.
Amplicon = namedtuple('Amplicon', ['sample', 'assay', 'locus', 'contig', 'strand',
                                   'start', 'end', 'size', 'f_mismatches', 'r_mismatches'])

COMPLEMENT = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')

FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.fas', '.seq', '.ffn']


class InSilicoPCR(object):

    def __init__(self, mismatches=2, three_prime=3, min_size=50, max_size=3000):
        """
        :param mismatches: type INT: Maximum number of mismatches tolerated per primer
        :param three_prime: type INT: Number of bases at the 3' end of each primer that must match exactly
        :param min_size: type INT: Minimum product size to report, primers included
        :param max_size: type INT: Maximum product size to report, primers included
        """
        self.mismatches = mismatches
        self.three_prime = three_prime
        self.min_size = min_size
        self.max_size = max_size

    @staticmethod
    def reverse_complement(seq):
        """
        :param seq: type STR: Nucleotide sequence, upper case
        :return: type STR: Reverse complement of seq
        """
        return seq.translate(COMPLEMENT)[::-1]

    @staticmethod
    def file_list(path):
        """
        Find the FASTA files to process
        :param path: type STR: Absolute path to a FASTA file or to a folder of FASTA files
        :return: type LIST: Sorted absolute paths of the FASTA files found
        """
        if os.path.isfile(path):
            return [os.path.abspath(path)]
        files = list()
        for entry in sorted(os.listdir(path)):
            name = entry[:-3] if entry.endswith('.gz') else entry
            if os.path.splitext(name)[1].lower() in FASTA_EXTENSIONS:
                files.append(os.path.join(os.path.abspath(path), entry))
        return files

    @staticmethod
    def strain_namer(fasta_file):
        """
        Derive a sample name from a file name, dropping the .gz and FASTA extensions
        :param fasta_file: type STR: Absolute path to a FASTA file
        :return: type STR: Sample name
        """
        name = os.path.basename(fasta_file)
        if name.endswith('.gz'):
            name = name[:-3]
        return os.path.splitext(name)[0]

    @staticmethod
    def read_fasta(fasta_file):
        """
        Read a (optionally gzipped) FASTA file into memory, upper cased
        :param fasta_file: type STR: Absolute path to a FASTA file
        :return: type LIST: (contig name, sequence) tuples
        """
        opener = gzip.open if fasta_file.endswith('.gz') else open
        with opener(fasta_file, 'rt') as handle:
            return [(record.id, str(record.seq).upper()) for record in SeqIO.parse(handle, 'fasta')]

    @staticmethod
    def seed_blocks(primer, max_mismatches):
        """
        Split a primer into max_mismatches + 1 contiguous blocks. A primer carrying no more than
        max_mismatches mismatches must match at least one block exactly, so the blocks can be used
        as exact-match seeds without missing any site.
        :param primer: type STR: Primer sequence
        :param max_mismatches: type INT: Maximum number of mismatches tolerated
        :return: type LIST: (offset within primer, block sequence) tuples
        """
        num_blocks = max_mismatches + 1
        length = len(primer)
        blocks = list()
        start = 0
        for index in range(num_blocks):
            stop = length * (index + 1) // num_blocks
            blocks.append((start, primer[start:stop]))
            start = stop
        return blocks

    def find_sites(self, seq, primer, anchor_three_prime):
        """
        Locate every position at which a primer matches the sequence within the allowed number of
        mismatches. Sites are seeded with exact block matches and then verified in full.
        :param seq: type STR: Nucleotide sequence to search, upper case
        :param primer: type STR: Primer sequence as searched, upper case
        :param anchor_three_prime: type BOOL: True when the 3' end of the primer is the right hand
        end of the searched sequence, False when the primer was reverse complemented and its 3' end
        is therefore the left hand end
        :return: type LIST: (start position, number of mismatches) tuples, 0-based, sorted
        """
        length = len(primer)
        if anchor_three_prime:
            exact = slice(length - self.three_prime, length)
        else:
            exact = slice(0, self.three_prime)
        exact_seq = primer[exact]
        # Seed on exact block matches to avoid scanning every position of the sequence
        candidates = set()
        for offset, block in self.seed_blocks(primer, self.mismatches):
            position = seq.find(block)
            while position != -1:
                start = position - offset
                if 0 <= start <= len(seq) - length:
                    candidates.add(start)
                position = seq.find(block, position + 1)
        # Verify each candidate site in full
        sites = list()
        for start in candidates:
            window = seq[start:start + length]
            if window[exact] != exact_seq:
                continue
            mismatches = 0
            for primer_base, seq_base in zip(primer, window):
                if primer_base != seq_base:
                    mismatches += 1
                    if mismatches > self.mismatches:
                        break
            if mismatches <= self.mismatches:
                sites.append((start, mismatches))
        return sorted(sites)

    def amplify_contig(self, seq, forward, reverse):
        """
        Find the products a primer pair yields on a single sequence, in both orientations
        :param seq: type STR: Nucleotide sequence to search, upper case
        :param forward: type STR: Forward primer, 5'->3'
        :param reverse: type STR: Reverse primer, 5'->3'
        :return: type LIST: (strand, start, end, forward mismatches, reverse mismatches) tuples,
        start and end 0-based and half open
        """
        products = dict()
        # A product is seen on the plus strand of the sequence as the forward primer followed by the
        # reverse complement of the reverse primer. When the locus sits the other way round in the
        # assembly, the same product is seen as the reverse primer followed by the reverse
        # complement of the forward primer.
        for strand, left, right in (('+', forward, self.reverse_complement(reverse)),
                                    ('-', reverse, self.reverse_complement(forward))):
            left_sites = self.find_sites(seq, left, anchor_three_prime=True)
            if not left_sites:
                continue
            right_sites = self.find_sites(seq, right, anchor_three_prime=False)
            if not right_sites:
                continue
            for left_start, left_mismatches in left_sites:
                # Pair each primer site with its nearest valid partner only. Pairing a site with
                # every downstream partner would invent long products across a repeated locus that
                # the short products outcompete in a real reaction.
                for right_start, right_mismatches in right_sites:
                    # Primers must not overlap
                    if right_start < left_start + len(left):
                        continue
                    end = right_start + len(right)
                    size = end - left_start
                    if size < self.min_size:
                        continue
                    # Sites are sorted, so every remaining partner is further away still
                    if size > self.max_size:
                        break
                    if strand == '+':
                        f_mismatches, r_mismatches = left_mismatches, right_mismatches
                    else:
                        f_mismatches, r_mismatches = right_mismatches, left_mismatches
                    # Keep the better supported product when a locus is found in both orientations
                    key = (left_start, end)
                    total = f_mismatches + r_mismatches
                    if key not in products or total < products[key][0]:
                        products[key] = (total, (strand, left_start, end, f_mismatches, r_mismatches))
                    break
        return [products[key][1] for key in sorted(products)]

    def amplify(self, fasta_file, assay, primers):
        """
        Run every primer pair of an assay against one FASTA file
        :param fasta_file: type STR: Absolute path to a FASTA file
        :param assay: type STR: Assay name
        :param primers: type DICT: Dictionary of locus name: {'F': forward primer, 'R': reverse primer}
        :return: type LIST: Amplicon named tuples
        """
        sample = self.strain_namer(fasta_file)
        amplicons = list()
        for contig, seq in self.read_fasta(fasta_file):
            for locus, pair in primers.items():
                for strand, start, end, f_mismatches, r_mismatches in self.amplify_contig(
                        seq=seq, forward=pair['F'], reverse=pair['R']):
                    amplicons.append(Amplicon(sample=sample,
                                              assay=assay,
                                              locus=locus,
                                              contig=contig,
                                              strand=strand,
                                              start=start + 1,
                                              end=end,
                                              size=end - start,
                                              f_mismatches=f_mismatches,
                                              r_mismatches=r_mismatches))
        return amplicons
