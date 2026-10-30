in_file = '/home/bioinfo/analyses/brucella_ladders/bruce_v2/consolidated_report/report_bruce.tsv'

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
            if species in filename and species not in ident and not any(ident == x for x in other_species):
                false_neg += 1
                notID += 1
            if species not in filename and species in ident:
                justPCR += 1
            if any(x in filename for x in other_species) and species in ident:
                false_pos += 1  # Assuming the preID is always the good one
            if species in filename and any(species == x for x in other_species):
                wrongPreID += 1  # Assuming the PCR ID is always the good one
            # elif species not in filename and any(ident == x for x in other_species):
            #     false_neg += 1

        print('Total ID (pre-ID or by PCR) {}: {}'.format(species, anyID))
        print('Genomes pre-ID as {}: {}'.format(species, preID))
        print('ID by PCR: {}/{} ({}%)'.format(pcrID, anyID, round(pcrID / anyID * 100, 1)))
        if pcrID > 0:
            print('Just ID by PCR: {}/{} ({}%)'.format(justPCR, pcrID, round(justPCR / pcrID * 100, 1)))
        else:
            print('Just ID by PCR: {}/{}'.format(justPCR, pcrID))
        print('Wrong pre-ID: {}/{} ({}%)'.format(wrongPreID, preID, round(wrongPreID / preID * 100, 1)))
        # print('True positives: {}/{} ({}%)'.format(true_pos, anyID, round(true_pos / anyID * 100, 1)))
        print('False positives: {}/{} ({}%)'.format(false_pos, preID, round(false_pos / preID * 100, 1)))
        print('False negatives (Pre-ID, but not ID by PCR, and not ID as other species): {}/{} ({}%)'.format(
            notID, anyID, round(notID / anyID * 100, 1)))
        # print('True negatives: {}/{} ({}%)'.format(true_neg, total, round(true_neg / total * 100, 1)))
        print('\n')
