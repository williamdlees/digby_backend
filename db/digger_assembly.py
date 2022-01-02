from receptor_utils import simple_bio_seq as simple
import os.path
import csv
from db.genomic_db import RefSeq, Feature, Sequence, Subject, SubjectSequence, Study
from db.genomic_db_functions import save_genomic_dataset_details, save_genomic_study, add_feature_to_ref, \
    save_genomic_sequence, save_genomic_ref_seq, feature_type, find_allele_by_seq, get_ref_set, find_or_assign_allele


def process_digger_record(session, study, assembly, dataset_dir, subject, annotation_file):
    print(f"Importing assembly {assembly.identifier} for subject {subject.identifier}")

    ref_seq = save_genomic_ref_seq(session, assembly.identifier, assembly.sequence, assembly.reference, assembly.chromosome, assembly.start, assembly.end)
    session.commit()

    # all records will be of the same sense. We'll use the sense when preparing the assembly file for gff

    sense = '+'

    with open(os.path.join(dataset_dir, annotation_file), 'r') as fi:
        reader = csv.DictReader(fi)
        feature_id = 1
        variants = {}

        for row in reader:
            if assembly.identifier in row['contig']:
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

                allele.features.append(region_feature)
                SubjectSequence(subject=subject, sequence=allele, haplotype='h1')
                # spruce up the ref if we have some info that can help
                if allele.imgt_name == '' and imgt_name != '':
                    allele.imgt_name = imgt_name

    session.commit()

    # Create assembly fasta for gff. The sequence name in the file must match the assembly id.

    ref_path = os.path.join(dataset_dir, f'{assembly.identifier}.fasta').replace(' ', '_')

    sequence = assembly.sequence if sense == '+' else simple.reverse_complement(assembly.sequence)

    with open(ref_path, 'w') as fo:
        fo.write(f'>{assembly.identifier}\n{assembly.sequence}\n')

    return
