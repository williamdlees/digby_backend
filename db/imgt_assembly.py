from lxml import html
import requests
from db.shared import delete_dependencies

from app import db
from db.feature_db import Species, RefSeq, Feature, Sequence, Sample, SampleSequence, Study
from db.save_genomic import save_genomic_dataset_details, save_genomic_study, save_genomic_sample, add_feature_to_ref, \
    save_genomic_sequence, save_genomic_ref_seq, feature_type


def process_imgt_assembly(sample_data):
    results = []
    needed_items = ['Assembly_id', 'Assembly_reference']

    valid_entry = True
    for needed_item in needed_items:
        if needed_item not in sample_data:
            results.append('%s not specified' % (needed_item))
            valid_entry = False

    if not valid_entry:
        return '\n'.join(results)

    # delete_dependencies(None)
    delete_dependencies(sample_data['Species'])

    results.append("Importing %s / %s" % (sample_data['Species'], sample_data['Sample']))

    page = requests.get(sample_data['URL'])
    tree = html.fromstring(page.content)

    seq_text = tree.xpath('//div[@class="sequence"]/pre')[0]

    sequence = ''
    for row in seq_text.text.split('\n'):
        if len(row) > 75:
            sequence += row[1:70].replace(' ', '')

    sp, data_set = save_genomic_dataset_details(sample_data['Locus'], sample_data['Dataset'], sample_data['Species'])
    ref_seq = save_genomic_ref_seq(sample_data['Locus'], sample_data['Assembly_id'], sp, sequence, sample_data['Assembly_reference'])
    db.session.commit()
    study = save_genomic_study(sample_data['Sample'], sample_data['Institute'], sample_data['Researcher'], sample_data['Reference'], sample_data['Contact'])

    sample = save_genomic_sample(sample_data['Sample'], sample_data['Type'], sample_data['Date'], study, sp.id, ref_seq.id, data_set.id, sample_data['URL'])

    features = tree.xpath('//div[@class="features"]/table')[0]
    rows = iter(features)

    state = None
    name = None
    gene_range = None
    strand = None
    parent_id = 0

    for row in rows:
        values = [col.text for col in row]

        if len(values) < 3:
            continue

        def get_range(s):
            gene_range = s.split('..')
            if len(gene_range) < 2:
                print('Invalid gene range found: %s' % s)
                return (('1', '1'), '+')

            for i in (0, 1):
                gene_range[i] = gene_range[i].replace('>', '').replace('<', '')

            strand = '+'

            if 'complement(' in gene_range[0]:
                gene_range[0] = gene_range[0].replace('complement(', '')
                gene_range[1] = gene_range[1].replace(')', '')
                strand = '-'

            try:
                if int(gene_range[1]) - int(gene_range[0]) > 10000000 or int(gene_range[0]) > int(gene_range[1]):
                    print('Invalid gene range found: %s' % s)
                    return ('1', '1'), '+'
            except:
                print('Invalid gene range found: %s' % s)
                return ('1', '1'), '+'

            return (gene_range, strand)

        if not state and values[0] in ['V-GENE', 'D-GENE', 'J-GENE']:
            gene_range, strand = get_range(values[2])
            state = values[0]

        elif state and not name:
            if values[1] == 'IMGT_allele':
                parent_id += 1
                name = values[2].split('*')[0]
                full_name = values[2]
                add_feature_to_ref(name, 'gene', gene_range[0], gene_range[1], strand, 'Name=%s;ID=%s' % (name, parent_id), parent_id, ref_seq)

                parent_id += 1
                add_feature_to_ref(name, 'mRNA', gene_range[0], gene_range[1], strand, 'Name=%s;ID=%s' % (name, parent_id), parent_id-1, ref_seq)

        elif state and name:
            if (state == 'V-GENE' and values[0] in ["5'UTR", 'L-PART1', 'V-INTRON', 'L-PART2', 'V-REGION', "3'UTR"]) \
                    or (state == 'D-GENE' and values[0] in ["5'UTR", 'D-REGION', "3'UTR"]) \
                    or (state == 'J-GENE' and values[0] in ["5'UTR", 'D-REGION', "3'UTR"]):
                gene_range, strand = get_range(values[2])

                if 'REGION' in values[0]:
                    seq_name = full_name
                    imgt_name = full_name
                else:
                    seq_name = name + '_' + values[0]
                    imgt_name = ''

                f = add_feature_to_ref(seq_name, feature_type[values[0]], gene_range[0], gene_range[1], strand, 'Name=%s;Parent=%s' % (name + '_' + values[0], parent_id), parent_id-1, ref_seq)
                s = save_genomic_sequence(seq_name, imgt_name, values[0], False, False, 'U', ref_seq.sequence[int(gene_range[0])-1:int(gene_range[1])], '', sp)

                s.features.append(f)
                SampleSequence(sample=sample, sequence=s, chromosome='h1,h2', chromo_count=2)

        if state and name and values[0] == "3'UTR":
            state = None
            name = None

    db.session.commit()
    return '\n'.join(results)
