from db.novel_alleles import find_novel_allele, assign_novel_gene
from db.shared import delete_dependencies
from Bio import SeqIO
import os.path
import csv
from app import db
from db.feature_db import Species, RefSeq, Feature, Sequence, Sample, SampleSequence, Study
from db.save_genomic import save_genomic_dataset_details, save_genomic_study, save_genomic_sample, add_feature_to_ref, \
    save_genomic_sequence, save_genomic_ref_seq, feature_type, find_allele_by_seq


def process_digger_record(sample_data, sample_dir, novels):
    results = []
    needed_items = ['Assembly_id', 'Assembly_reference', 'Annotation_file', 'Assembly_file', 'Chromosome', 'Start_CoOrd', 'End_CoOrd', 'Sample_description', 'Study_description']

    valid_entry = True
    for needed_item in needed_items:
        if needed_item not in sample_data:
            results.append('%s not specified' % (needed_item))
            valid_entry = False

    if not valid_entry:
        return '\n'.join(results)

    results.append("Importing %s / %s" % (sample_data['Species'], sample_data['Sample']))

    if '>' not in sample_data['Assembly_file']:
        sequence = SeqIO.read(os.path.join(sample_dir, sample_data['Assembly_file']), 'fasta')
    else:
        sequence = None
        seq_file, seq_name = os.path.join(sample_dir, sample_data['Assembly_file']).split('>')
        sequences = SeqIO.parse(seq_file, 'fasta')
        for record in sequences:
            if seq_name in record.description:
                sequence = record
                break
        if sequence is None:
            results.append('Sequence %s not found in file %s' % (seq_name, seq_file))
            return '\n'.join(results)

    sp, data_set = save_genomic_dataset_details(sample_data['Locus'], sample_data['Dataset'], sample_data['Species'])
    ref_seq = save_genomic_ref_seq(sample_data['Locus'], sample_data['Assembly_id'], sp, sequence, sample_data['Assembly_reference'],
                                   sample_data['Chromosome'], sample_data['Start_CoOrd'], sample_data['End_CoOrd'])
    db.session.commit()
    study = save_genomic_study(sample_data['Sample'], sample_data['Institute'], sample_data['Researcher'], sample_data['Reference'], sample_data['Contact'], sample_data['Study_description'])

    contig_filter = None
    annotation_filename = sample_data['Annotation_file']

    if '>' in annotation_filename:
        annotation_filename, contig_filter = annotation_filename.split('>')

    report_link = '/'.join(['study_data', 'Genomic', sample_data['Species'].replace(' ', '_'), sample_data['Dataset'].replace(' ', '_'), annotation_filename])
    sample = save_genomic_sample(sample_data['Sample'], sample_data['Type'], sample_data['Date'], study, sp.id, ref_seq.id, data_set.id, report_link, sample_data['Sample_description'])

    # all records will be of the same sense. We'll use the sense when preparing the assembly file for gff

    sense = '+'

    with open(os.path.join(sample_dir, annotation_filename), 'r') as fi:
        reader = csv.DictReader(fi)
        feature_id = 1
        variants = {}

        for row in reader:
            if contig_filter is None or contig_filter in row['contig']:
                sense = row['sense']

                # swap forward and reverse co-ords as analysis is on anti-sense

                row['start'] = row['start_rev']
                row['end'] = row['end_rev']
                for k in row.keys():
                    if '_start_rev' in k:
                        pos_name = k.replace('_start_rev', '')
                        row[pos_name + '_start'] = row[pos_name + '_start_rev']
                        row[pos_name + '_end'] = row[pos_name + '_end_rev']

                imgt_name = ''
                if len(row['imgt_score']) > 0 and float(row['imgt_score']) >= 99.9:
                    imgt_name = row['imgt_match']

                allele = find_allele_by_seq(row['seq'], sp.id)

                if allele is not None:
                    name = allele.name
                elif len(row['cirelli_score']) > 0 and float(row['cirelli_score']) > 99.99:
                    name = row['cirelli_match']
                elif len(row['cirelli_match']) > 0:
                    allele_info = find_novel_allele(novels, sample_data['Species'], sample_data['Dataset'], row['seq'])

                    if allele_info is None:
                        # we know we have a cirelli name: use it for family and locus
                        stem = row['cirelli_match'].split('_')[1]
                        stem = stem.split('.')[0]
                        family = stem[-1:]
                        prefix = 'VDJbase_' + stem[:-1]
                        allele_info = assign_novel_gene(novels, sample_data['Species'], sample_data['Dataset'], row['seq'], prefix, family, 'Digger')

                    name = allele_info['name']
                else:
                    continue    # poor sequence, probably truncated alignment

                add_feature_to_ref(name, 'gene', row['start'], row['end'], row['sense'], 'Name=%s;ID=%s' % (name, feature_id), feature_id, ref_seq)
                parent_id = feature_id
                feature_id += 1

                add_feature_to_ref(name, 'mRNA', row['start'], row['end'], row['sense'], 'Name=%s;ID=%s' % (name, parent_id), parent_id-1, ref_seq)

                if 'IGHV' in name:
                    allele_type = 'V-REGION'
                    region_feature = add_feature_to_ref(name, 'CDS', row['start'], row['end'], row['sense'], 'Name=%s_V_REGION;ID=%s' % (name, feature_id), parent_id, ref_seq)
                    feature_id += 1
                    if row['3_rss_start']:
                        add_feature_to_ref(name, 'CDS', row['3_rss_start'], row['3_rss_end'], row['sense'], 'Name=%s_V-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1
                    if row['l_part1_start']:
                        add_feature_to_ref(name, 'CDS', row['l_part1_start'], row['l_part1_end'], row['sense'], 'Name=%s_L-PART1;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1
                    if row['l_part2_start']:
                        add_feature_to_ref(name, 'CDS', row['l_part2_start'], row['l_part2_end'], row['sense'], 'Name=%s_L-PART2;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1

                elif 'IGHD' in name:
                    allele_type = 'D-REGION'
                    region_feature = add_feature_to_ref(name, 'CDS', row['start'], row['end'], row['sense'], 'Name=%s_D_REGION;ID=%s' % (name, feature_id), parent_id, ref_seq)
                    feature_id += 1
                    if row['3_rss_start']:
                        add_feature_to_ref(name, 'CDS', row['3_rss_start'], row['3_rss_end'], row['sense'], 'Name=%s_3D-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1
                    if row['5_rss_start']:
                        add_feature_to_ref(name, 'CDS', row['5_rss_start'], row['5_rss_end'], row['sense'], 'Name=%s_5D-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1

                elif 'IGHJ' in name:
                    region_feature = add_feature_to_ref(name, 'CDS', row['start'], row['end'], row['sense'], 'Name=%s_J_REGION;ID=%s' % (name, feature_id), parent_id, ref_seq)
                    feature_id += 1
                    if row['5_rss_start']:
                        add_feature_to_ref(name, 'CDS', row['5_rss_start'], row['5_rss_end'], row['sense'], 'Name=%s_J-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        allele_type = 'J-REGION'

                if allele is None:
                    functional = 'U'
                    if row['functional'] == 'functional':
                        functional = 'F'
                    elif row['functional'] == 'pseudo':
                        functional = 'P'
                    elif row['functional'] == 'ORF':
                        functional = 'O'
                    gs = save_genomic_sequence(name, imgt_name, allele_type, True, False, functional, row['seq'], row['seq_gapped'], sp)
                    gs.features.append(region_feature)
                    SampleSequence(sample=sample, sequence=gs, chromosome='h1', chromo_count=1)
                else:
                    allele.features.append(region_feature)
                    SampleSequence(sample=sample, sequence=allele, chromosome='h1', chromo_count=1)
                    # spruce up the ref if we have some info that can help
                    if allele.gapped_sequence == '':
                        allele.gapped_sequence = row['seq_gapped'].lower()
                    if allele.imgt_name == '' and imgt_name != '':
                        allele.imgt_name = imgt_name

    db.session.commit()

    # Create assembly fasta for gff. The sequence name in the file must match the assembly id.

    ref_path = os.path.join(sample_dir, '%s_%s.fasta' % (sample_data['Species'], sample_data['Assembly_id'])).replace(' ', '_')

    if sense == '-':
        sequence = sequence.reverse_complement()

    with open(ref_path, 'w') as fo:
        fo.write('>%s\n' % sample_data['Assembly_id'])
        fo.write(str(sequence.seq))


    return '\n'.join(results)
