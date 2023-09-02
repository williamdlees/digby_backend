from receptor_utils import simple_bio_seq as simple
import os.path
import shutil
import csv

from db.genomic_db import RefSeq, Feature, Sequence, Subject, SampleSequence, Study
from db.genomic_db_functions import save_genomic_dataset_details, save_genomic_study, add_feature_to_ref, \
    save_genomic_sequence, save_genomic_ref_seq, find_allele_by_seq, get_ref_set, find_or_assign_allele, link_sequence_to_feature, update_sample_sequence_link


def process_digger_record(session, species, assembly, dataset_dir, subject, annotation_file, reference_features):
    print(f"Importing assembly {assembly.identifier} for subject {subject.identifier}")

    ref_seq = save_genomic_ref_seq(session, assembly.identifier, assembly.sequence, assembly.reference, assembly.chromosome, assembly.start, assembly.end)

    shutil.copy(annotation_file, os.path.join(dataset_dir, 'samples'))

    # all records will be of the same sense. We'll use the sense when preparing the assembly file for gff

    with open(os.path.join(dataset_dir, annotation_file), 'r') as fi:
        reader = csv.DictReader(fi)
        feature_id = 1
        variants = {}

        for row in reader:
            if assembly.identifier in row['contig']:
                # all records will be of the same sense. We'll use the sense when preparing the assembly file for gff
                sense = row['sense']

                # swap forward and reverse co-ords as analysis is on anti-sense
                if sense == '-':
                    row['end'] = row['start_rev']
                    row['start'] = row['end_rev']
                    for k in row.keys():
                        if '_start_rev' in k:
                            pos_name = k.replace('_start_rev', '')
                            row[pos_name + '_end'] = row[pos_name + '_start_rev']
                            row[pos_name + '_start'] = row[pos_name + '_end_rev']

                imgt_name = ''
                if len(row['imgt_score']) > 0 and float(row['imgt_score']) >= 99.9:
                    imgt_name = row['imgt_match']

                allele = find_allele_by_seq(session, row['seq'])

                if allele is None:
                    if len(row['cirelli_match']) > 0:
                        functional = 'U'
                        if row['functional'].lower() == 'functional':
                            functional = 'F'
                        elif row['functional'] == 'pseudo':
                            functional = 'P'
                        elif row['functional'] == 'ORF':
                            functional = 'O'
                        allele = find_or_assign_allele(session, row['seq'], 'V' in row['cirelli_match'], functional)
                    else:
                        continue    # poor sequence, probably truncated alignment

                name = allele.name

                def add_feature(name, feature_type, feature_seq, feature, start, end, strand, attribute, parent_id, ref):
                    start = int(start)
                    end = int(end)
                    add_feature_to_ref(name, 'gene', feature_type, feature_seq, '', 'gene', start, end, strand, attribute, parent_id, ref)
                    add_feature_to_ref(name, 'gene', feature_type, feature_seq, '', 'mRNA', start, end, strand, attribute, parent_id, ref)
                    return add_feature_to_ref(name, 'allele', feature_type, feature_seq, '', feature, start, end, strand, attribute, parent_id, ref)

                if 'IGHV' in name:
                    region_feature = add_feature(name, 'V-REGION', row['seq'], 'CDS', row['start'], row['end'], row['sense'], 'Name=%s_V_REGION;ID=%s' % (name, feature_id), feature_id, ref_seq)
                    parent_id = feature_id
                    feature_id += 1
                    if row['3_rss_start']:
                        add_feature(name, 'V-RS', row['seq'], 'CDS', row['3_rss_start'], row['3_rss_end'], row['sense'], 'Name=%s_V-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1
                    if row['l_part1_start']:
                        add_feature(name, 'L-PART1', row['seq'], 'CDS', row['l_part1_start'], row['l_part1_end'], row['sense'], 'Name=%s_L-PART1;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1
                    if row['l_part2_start']:
                        add_feature(name, 'L-PART2', row['seq'], 'CDS', row['l_part2_start'], row['l_part2_end'], row['sense'], 'Name=%s_L-PART2;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1

                elif 'IGHD' in name:
                    region_feature = add_feature(name, 'D-REGION', row['seq'], 'CDS', row['start'], row['end'], row['sense'], 'Name=%s_D_REGION;ID=%s' % (name, feature_id), feature_id, ref_seq)
                    parent_id = feature_id
                    feature_id += 1
                    if row['3_rss_start']:
                        add_feature(name, '3D-RS', row['seq'], 'CDS', row['3_rss_start'], row['3_rss_end'], row['sense'], 'Name=%s_3D-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1
                    if row['5_rss_start']:
                        add_feature(name, '5D-RS', row['seq'], 'CDS', row['5_rss_start'], row['5_rss_end'], row['sense'], 'Name=%s_5D-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)
                        feature_id += 1

                elif 'IGHJ' in name:
                    region_feature = add_feature(name, 'J-REGION', row['seq'], 'CDS', row['start'], row['end'], row['sense'], 'Name=%s_J_REGION;ID=%s' % (name, feature_id), feature_id, ref_seq)
                    feature_id += 1
                    parent_id = feature_id
                    if row['5_rss_start']:
                        add_feature(name, 'J-RS', row['seq'], 'CDS', row['5_rss_start'], row['5_rss_end'], row['sense'], 'Name=%s_J-RS;ID=%s' % (name, feature_id), parent_id, ref_seq)

                else:
                    breakpoint()

                allele.features.append(region_feature)
                # spruce up the ref if we have some info that can help
                if allele.imgt_name == '' and imgt_name != '':
                    allele.imgt_name = imgt_name

                update_subject_sequence_link(session, 0, subject, allele)


    # Create assembly fasta for gff. The sequence name in the file must match the assembly id.

    ref_path = os.path.join(dataset_dir, 'samples', f"{species.replace(' ', '_')}_{assembly.identifier}.fasta")

    sequence = assembly.sequence if sense == '+' else simple.reverse_complement(assembly.sequence)

    with open(ref_path, 'w') as fo:
        fo.write(f'>{assembly.identifier}\n{sequence}\n')

    return
