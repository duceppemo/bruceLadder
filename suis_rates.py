import sys

# Path to the two column report_suis.tsv written by bruce_ladder.py, as the first argument
in_file = sys.argv[1] if len(sys.argv) > 1 else 'report_suis.tsv'

biovar_list = ['bv01', 'bv02', 'bv03', 'bv04', 'bv05']

for biovar in biovar_list:
    if biovar == 'bv04':
        t=1
    bv = biovar[-1]
    preID = 0
    pos = 0
    pcrID = 0
    anyID = 0
    justPCR = 0
    false_pos = 0
    false_neg = 0
    with open(in_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line == '':
                continue

            filename, ident = line.split('\t')

            if 'suis' in filename and 'bv-' + bv in filename:
                preID += 1
            if biovar in ident:
                pcrID += 1
            if ('suis' in filename and 'bv-' + bv in filename) or biovar in ident:
                anyID += 1

            if 'suis' not in filename and biovar in ident:
                justPCR += 1

            if 'suis' in filename and 'bv-' + bv in filename and biovar in ident:
                pos += 1
            if 'suis' in filename and 'bv-' + bv in filename and biovar not in ident:
                false_neg += 1
            if 'suis' not in filename and biovar in ident:
                false_pos += 1

        print('Total ID (pre-ID or by PCR) {}: {}'.format(biovar, anyID))
        print('Genomes pre-ID as suis {}: {}'.format(biovar, preID))
        if anyID > 0:
            print('ID by PCR: {}/{} ({}%)'.format(pcrID, anyID, round(pcrID / anyID * 100, 1)))
        else:
            print('ID by PCR: {}/{}'.format(pcrID, anyID))
        if pcrID > 0:
            print('Just ID by PCR: {}/{} ({}%)'.format(justPCR, pcrID, round(justPCR / pcrID * 100, 1)))
        else:
            print('Just ID by PCR: {}/{}'.format(justPCR, pcrID))

        if preID > 0:
            # print('Correct biovar ID: {}/{} ({}%)'.format(pos, pcrID, round(pos / pcrID * 100, 1)))
            print('False positives: {}/{} ({}%)'.format(false_pos, preID, round(false_pos / preID * 100, 1)))
            print('False negatives: {}/{} ({}%)'.format(false_neg, anyID, round(false_neg / anyID * 100, 1)))
        print('\n')
