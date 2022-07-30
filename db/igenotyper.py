import re
from hashlib import sha256

from receptor_utils import simple_bio_seq as simple
import os.path
import shutil
import csv

from receptor_utils.novel_allele_name import name_novel
from receptor_utils.simple_bio_seq import write_csv
from sqlalchemy import and_

from db.genomic_db import RefSeq, Feature, Sequence, Subject, SubjectSequence, Study
from db.genomic_db_functions import rationalise_name, save_novel_allele, find_gene_allele_by_seq, add_feature_to_ref, find_feature_by_name, \
    find_feature_by_sequence, find_sequence_by_sequence, save_genomic_sequence, link_sequence_to_feature, update_subject_sequence_link, find_allele_by_name


class GeneParsingException(Exception):
    pass

annotation_records = {}

def read_csv(filename):
    with open(filename, 'r') as fi:
        return list(csv.DictReader(fi))

def process_igenotyper_record(session, species, dataset_dir, subject, annotation_file, reference_features):
    global annotation_records

    print(f"Importing subject {subject.identifier}")

    if annotation_file not in annotation_records:
        annotation_records[annotation_file] = read_csv(annotation_file)

    rows = [x for x in annotation_records[annotation_file] if x['sample_name'] == subject.name_in_study]

    # If the annotation file contains records for multiple subjects, split into multiple files

    if len(rows) < len(annotation_records[annotation_file]):
        subj_af_name, subj_af_ext = os.path.splitext(os.path.basename(annotation_file))
        subj_af_name = f"{subj_af_name}_{subject.identifier}{subj_af_ext}"
        write_csv(os.path.join(dataset_dir, 'samples', subj_af_name), rows)
        subject.annotation_path = subject.annotation_path.replace(os.path.basename(annotation_file), subj_af_name)
    elif not os.path.isfile(os.path.join(dataset_dir, 'samples')):
        shutil.copy(annotation_file, os.path.join(dataset_dir, 'samples'))

    sense = '+'     # by + sense we mean 5' to 3'
    feature_id = 1

    for row in rows:
        gene_type = row['genotyper_gene'][3]

        if gene_type in ['V', 'D', 'J']:
            seq = find_allele_by_name(session, row['vdjbase_allele'])

            if not seq:
                if gene_type == 'V':
                    gapped_seq = row['V-REGION-GAPPED']
                else:
                    gapped_seq = ''

                seq = save_novel_allele(session, row['genotyper_gene'], row['vdjbase_allele'], row['notes'].replace('\\n', '\r\n'), row[f'{gene_type}-REGION'], gapped_seq)

            update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, seq)
            feature = find_feature_by_name(session, f'{gene_type}-REGION', seq.name, subject.ref_seq)

            if feature and seq.sequence.replace('.', '') != feature.feature_seq.replace('.', ''):
                print(f'Error: feature {feature.name} sequence does not match that of sequence {seq.name} in subject {subject.identifier}')

            if gene_type == 'V':
                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[subject.ref_seq.name][seq.gene.name]['exon_2']['start'] + 11
                    end = reference_features[subject.ref_seq.name][seq.gene.name]['exon_2']['end']

                    # if there are gaps in the alignment against the reference assembly, carry these into the feature

                    if '.' in row['V-REGION']:
                        feature_seq = row['V-REGION']
                        seq.sequence = row['V-REGION']
                    else:
                        feature_seq = seq.sequence

                    feature = add_feature_to_ref(seq.name, 'allele', 'V-REGION', feature_seq, row['V-REGION-CIGAR'], 'CDS', start, end, '+',
                                                 f"Name={seq.name}_V-REGION;ID={feature_id}", feature_id, subject.ref_seq)

                link_sequence_to_feature(seq, feature)
                add_feature('V-NONAMER', 'nonamer', reference_features, row, seq, session, subject)
                add_feature('V-SPACER', 'spacer', reference_features, row, seq, session, subject)
                add_feature('V-HEPTAMER', 'heptamer', reference_features, row, seq, session, subject)
                add_feature('L-PART2', 'exon_2', reference_features, row, seq, session, subject)
                add_feature('V-INTRON', 'gencode_intron', reference_features, row, seq, session, subject)
                add_feature('L-PART1', 'exon_1', reference_features, row, seq, session, subject)

            elif gene_type == 'D':
                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[subject.ref_seq.name][seq.gene.name]['exon_1']['start']
                    end = reference_features[subject.ref_seq.name][seq.gene.name]['exon_1']['end']
                    feature = add_feature_to_ref(seq.name, 'allele', 'D-REGION', seq.sequence, row['D-REGION-CIGAR'], 'CDS', start, end, '+',
                                                 f"Name={seq.name}_D-REGION;ID={feature_id}", feature_id, subject.ref_seq)

                link_sequence_to_feature(seq, feature)
                add_feature('D-3_NONAMER', '3_nonamer', reference_features, row, seq, session, subject)
                add_feature('D-3_SPACER', '3_spacer', reference_features, row, seq, session, subject)
                add_feature('D-3_HEPTAMER', '3_heptamer', reference_features, row, seq, session, subject)
                add_feature('D-5_NONAMER', '5_nonamer', reference_features, row, seq, session, subject)
                add_feature('D-5_SPACER', '5_spacer', reference_features, row, seq, session, subject)
                add_feature('D-5_HEPTAMER', '5_heptamer', reference_features, row, seq, session, subject)

            elif gene_type == 'J':
                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[subject.ref_seq.name][seq.gene.name]['exon_1']['start']
                    end = reference_features[subject.ref_seq.name][seq.gene.name]['exon_1']['end']
                    feature = add_feature_to_ref(seq.name, 'allele', 'J-REGION', seq.sequence, row['J-REGION-CIGAR'], 'CDS', start, end, '+',
                                                 f"Name={seq.name}_J-REGION;ID={feature_id}", feature_id, subject.ref_seq)

                link_sequence_to_feature(seq, feature)
                add_feature('J-NONAMER', 'nonamer', reference_features, row, seq, session, subject)
                add_feature('J-SPACER', 'spacer', reference_features, row, seq, session, subject)
                add_feature('J-HEPTAMER', 'heptamer', reference_features, row, seq, session, subject)

    session.commit()



def add_feature(feature, bed_name, reference_features, row, seq, session, subject):
    if not row[feature]:
        return

    feature_seq = find_sequence_by_sequence(session, feature, seq.gene.name, row[feature])

    if not feature_seq:
        feature_name = f"{seq.gene.name}*{sha256(row[feature].encode('utf-8')).hexdigest()[-4:]}"
        feature_seq = save_genomic_sequence(session, feature_name, seq.gene.name, feature, True, False, '', row[feature], '')

    update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, feature_seq)

    feature_rec = find_feature_by_name(session, feature, feature_seq.name, subject.ref_seq)

    if not feature_rec:
        feature_id = session.query(Feature).count()

        start = reference_features[subject.ref_seq.name][seq.gene.name][bed_name]['start']

        if feature != 'L-PART2':
            end = reference_features[subject.ref_seq.name][seq.gene.name][bed_name]['end']
        else:
            end = reference_features[subject.ref_seq.name][seq.gene.name][bed_name]['start'] + 11

        feature_rec = add_feature_to_ref(feature_seq.name, 'allele', feature, feature_seq.sequence, '', 'UTR', start, end, '+',
                                     f"Name={feature_seq.name};ID={feature_id}", feature_id, subject.ref_seq)

    link_sequence_to_feature(feature_seq, feature_rec)


def add_gene_level_subfeature(feature, imgt_feature_name, name_prefix, feature_id, parent_id, ref):
    if imgt_feature_name == 'exon_2':
        name = f"IGHVRegion{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
        add_feature_to_ref(name, 'gene', 'V-REGION', feature['ref_seq'], '', 'CDS', feature['start']+11, feature['end'], '+',
                           f"Name={feature['gene']}_V-REGION;ID={feature_id}", parent_id, ref)

        name = f"IGHVLpart2{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
        add_feature_to_ref(name, 'gene', 'L_PART-2', feature['ref_seq'], '', 'CDS', feature['start'], feature['start']+11, '+',
                           f"Name={feature['gene']}_L_PART2;ID={feature_id}", parent_id, ref)
    else:
        name = f"{name_prefix}{feature['gene'][3:]}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
        add_feature_to_ref(name, 'gene', imgt_feature_name, feature['ref_seq'], '', 'CDS', feature['start'], feature['end'], '+',
                           f"Name={feature['gene']}_{imgt_feature_name};ID={feature_id}", parent_id, ref)



def add_gene_level_features(session, ref, reference_features):
    feature_id = 1
    for gene, features in reference_features[ref.name].items():
        parent_id = feature_id
        for feature_type, feature in features.items():
            if 'V' in feature['gene']:
                if feature_type == 'gene':
                    add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                elif feature_type == 'nonamer':
                    add_gene_level_subfeature(feature, 'V-NONAMER', 'IGHVNona', feature_id, parent_id, ref)
                elif feature_type == 'spacer':
                    add_gene_level_subfeature(feature, 'V-SPACER', 'IGHVSpacer', feature_id, parent_id, ref)
                elif feature_type == 'heptamer':
                    add_gene_level_subfeature(feature, 'V-HEPTAMER', 'IGHVHepta', feature_id, parent_id, ref)
                    # TODO - split bed exon_2 into L_PART2, V-REGION
                elif feature_type == 'exon_2':
                    add_gene_level_subfeature(feature, 'exon_2', 'IGHVRegion', feature_id, parent_id, ref)
                elif feature_type == 'gencode_intron':
                    add_gene_level_subfeature(feature, 'V-INTRON', 'IGHVIntron', feature_id, parent_id, ref)
                elif feature_type == 'exon_1':
                    add_gene_level_subfeature(feature, 'L_PART-1', 'IGHVLP1', feature_id, parent_id, ref)

            elif 'D' in feature['gene']:
                if feature_type == 'gene':
                    add_feature_to_ref(feature['gene'], 'gene', 'D-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'D-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                elif feature_type == '3_nonamer':
                    add_gene_level_subfeature(feature, 'D-3-NONAMER', 'IGHD3Nona', feature_id, parent_id, ref)
                elif feature_type == '3_spacer':
                    add_gene_level_subfeature(feature, 'D-3-SPACER', 'IGHD3Spacer', feature_id, parent_id, ref)
                elif feature_type == '3_heptamer':
                    add_gene_level_subfeature(feature, 'D-3-HEPTAMER', 'IGHD3Hepta', feature_id, parent_id, ref)
                elif feature_type == '5_nonamer':
                    add_gene_level_subfeature(feature, 'D-5-NONAMER', 'IGHD5Nona', feature_id, parent_id, ref)
                elif feature_type == '5_spacer':
                    add_gene_level_subfeature(feature, 'D-5-SPACER', 'IGHD5Spacer', feature_id, parent_id, ref)
                elif feature_type == '5_heptamer':
                    add_gene_level_subfeature(feature, 'D-5-HEPTAMER', 'IGHD5Hepta', feature_id, parent_id, ref)
                elif feature_type == 'exon_1':
                    add_gene_level_subfeature(feature, 'D-REGION', 'IGHDRegion', feature_id, parent_id, ref)

            elif 'J' in feature['gene']:
                if feature_type == 'gene':
                    add_feature_to_ref(feature['gene'], 'gene', 'J-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'J-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                elif feature_type == 'nonamer':
                    add_gene_level_subfeature(feature, 'J-NONAMER', 'IGHJNona', feature_id, parent_id, ref)
                elif feature_type == 'spacer':
                    add_gene_level_subfeature(feature, 'J-SPACER', 'IGHJSpacer', feature_id, parent_id, ref)
                elif feature_type == 'heptamer':
                    add_gene_level_subfeature(feature, 'J-HEPTAMER', 'IGHJHepta', feature_id, parent_id, ref)
                elif feature_type == 'exon_1':
                    add_gene_level_subfeature(feature, 'J-REGION', 'IGHJRegion', feature_id, parent_id, ref)

            feature_id += 1



