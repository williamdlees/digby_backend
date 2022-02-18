# Read a bed file

import csv

def read_bed_file(infile):
    r_headers = [
        'ref_name',
        'start',
        'end',
        'gene',
        'feature',
        'ref_seq',
    ]

    with open(infile, 'r') as fi:
        reader = csv.DictReader(fi, delimiter='\t', fieldnames=r_headers)
        recs = {}
        for row in reader:
            if row['feature'] is None and row['ref_seq'] is None:
                row['feature'] = infile.split('.')[0]
                if row['feature'] == 'genes':
                    row['feature'] = 'gene'
                row['ref_seq'] = ''
            if row['ref_name'] not in recs:
                recs[row['ref_name']] = []
            recs[row['ref_name']].append(row)
        return recs


# Read multiple bed files into a single set organised by ref file, gene and feature
# convert co-ords to + sense if needed


def read_bed_files(bed_files, sense, ref_seq_len):
    features = {}

    for infile in bed_files:
        bed_data = read_bed_file(infile)

        for ref_file, recs in bed_data.items():
            if ref_file not in features:
                features[ref_file] = {}
            for rec in recs:
                if rec['start'] and rec['end']:
                    rec['start'] = int(rec['start'])
                    rec['end'] = int(rec['end'])

                    if sense == '-':
                        rec['start'], rec['end'] = ref_seq_len - rec['end'] + 1, ref_seq_len - rec['start']

                    if rec['gene'] not in features[ref_file]:
                        features[ref_file][rec['gene']] = {}
                    if rec['feature'] not in features[ref_file][rec['gene']]:
                        features[ref_file][rec['gene']][rec['feature']] = {}
                    features[ref_file][rec['gene']][rec['feature']] = rec

    return features
