import sys

# Path to the two column report_bruce.tsv written by bruce_ladder.py, as the first argument
in_file = sys.argv[1] if len(sys.argv) > 1 else 'report_bruce.tsv'


def fraction(numerator, denominator):
    # A species that is absent from the report leaves its totals at zero, so the percentage is
    # dropped rather than divided
    if denominator == 0:
        return '{}/{}'.format(numerator, denominator)
    return '{}/{} ({}%)'.format(numerator, denominator, round(numerator / denominator * 100, 1))


def called_as(ident, names):
    # The species_list entries are substrings of the calls bruce_ladder.py writes, so a name has
    # to be looked for inside the call: 'abort' inside 'abortus', 'suis' inside
    # 'suis_bv2/suis_bv1_3_4'. Comparing the two for equality misses those.
    return any(name in ident for name in names)


species_list = ['abort', 'suis', 'ovis', 'canis', 'ceti', 'pedialis', 'melitensis', 'neotomae']

for species in species_list:
    preID = 0
    pcrID = 0
    anyID = 0
    notID = 0
    justPCR = 0
    wrongPreID = 0
    true_pos = 0
    false_pos = 0
    true_neg = 0
    false_neg = 0
    total = 0

    other_species = species_list.copy()
    other_species.remove(species)
    with open(in_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line == '':
                continue

            filename, ident = line.split('\t')

            total += 1
            if species in filename:
                preID += 1
            if species in ident:
                pcrID += 1
            if species in filename or species in ident:
                anyID += 1

            if species in filename and species in ident:
                true_pos += 1
            if species in filename and species not in ident and not called_as(ident, other_species):
                false_neg += 1
                notID += 1
            if species not in filename and species in ident:
                justPCR += 1
            if any(x in filename for x in other_species) and species in ident:
                false_pos += 1  # Assuming the preID is always the good one
            if species in filename and species not in ident and called_as(ident, other_species):
                wrongPreID += 1  # Assuming the PCR ID is always the good one
            # elif species not in filename and any(ident == x for x in other_species):
            #     false_neg += 1

        print('Total ID (pre-ID or by PCR) {}: {}'.format(species, anyID))
        print('Genomes pre-ID as {}: {}'.format(species, preID))
        print('ID by PCR: {}'.format(fraction(pcrID, anyID)))
        print('Just ID by PCR: {}'.format(fraction(justPCR, pcrID)))
        print('Wrong pre-ID: {}'.format(fraction(wrongPreID, preID)))
        # print('True positives: {}'.format(fraction(true_pos, anyID)))
        print('False positives: {}'.format(fraction(false_pos, preID)))
        print('False negatives (Pre-ID, but not ID by PCR, and not ID as other species): {}'.format(
            fraction(notID, anyID)))
        # print('True negatives: {}'.format(fraction(true_neg, total)))
        print('\n')
