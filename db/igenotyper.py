import re
from hashlib import sha256

from receptor_utils import simple_bio_seq as simple
import os.path
import shutil
import csv

from receptor_utils.novel_allele_name import name_novel
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

    if not os.path.isfile(os.path.join(dataset_dir, 'samples')):
        shutil.copy(annotation_file, os.path.join(dataset_dir, 'samples'))

    if annotation_file not in annotation_records:
        annotation_records[annotation_file] = read_csv(annotation_file)

    rows = [x for x in annotation_records[annotation_file] if x['sample_name'] == subject.name_in_study]

    sense = '+'     # by + sense we mean 5' to 3'
    feature_id = 1

    for row in rows:
        chain = row['genotyper_gene'][3]

        if chain == 'V':
            seq = find_allele_by_name(session, row['vdjbase_allele'])

            if not seq:
                seq = save_novel_allele(session, row['genotyper_gene'], row['vdjbase_allele'], row['notes'].replace('\\n', '\r\n'), row['V-REGION'], row['V-REGION-GAPPED'])

            update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, seq)
            feature = find_feature_by_name(session, 'V-REGION', seq.name, subject.ref_seq)

            if feature and seq.sequence != feature.feature_seq:
                print(f'Error: feature {feature.name} sequence does not match that of sequence {seq.name} in subject {subject.identifier}')

            if not feature:
                feature_id = session.query(Feature).count()
                start = reference_features[subject.ref_seq.name][seq.gene]['exon_2']['start'] + 11
                end = reference_features[subject.ref_seq.name][seq.gene]['exon_2']['end']
                feature = add_feature_to_ref(seq.name, 'allele', 'V-REGION', seq.sequence, 'CDS', start, end, '+',
                               f"Name={seq.name}_V-REGION;ID={feature_id}", feature_id, subject.ref_seq)

            link_sequence_to_feature(seq, feature)

            # TODO - use IMGT names in bed files and get rid of one of these arguments

            add_feature('V-NONAMER', 'nonamer', reference_features, row, seq, session, subject)
            add_feature('V-SPACER', 'spacer', reference_features, row, seq, session, subject)
            add_feature('V-HEPTAMER', 'heptamer', reference_features, row, seq, session, subject)
            add_feature('L-PART2', 'exon_2', reference_features, row, seq, session, subject)
            add_feature('V-INTRON', 'gencode_intron', reference_features, row, seq, session, subject)
            add_feature('L-PART1', 'exon_1', reference_features, row, seq, session, subject)


def add_feature(feature, bed_name, reference_features, row, seq, session, subject):
    if not row[feature]:
        return

    feature_seq = find_sequence_by_sequence(session, feature, seq.gene, row[feature])

    if not feature_seq:
        feature_name = f"{seq.gene}*{sha256(row[feature].encode('utf-8')).hexdigest()[-4:]}"
        feature_seq = save_genomic_sequence(session, feature_name, seq.gene, feature, True, False, '', row[feature], '')

    update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, feature_seq)

    feature_rec = find_feature_by_name(session, feature, feature_seq.name, subject.ref_seq)

    if not feature_rec:
        feature_id = session.query(Feature).count()
        start = reference_features[subject.ref_seq.name][seq.gene][bed_name]['start']

        if feature != 'L-PART2':
            end = reference_features[subject.ref_seq.name][seq.gene][bed_name]['end']
        else:
            end = reference_features[subject.ref_seq.name][seq.gene][bed_name]['start'] + 11

        feature_rec = add_feature_to_ref(feature_seq.name, 'allele', feature, feature_seq.sequence, 'UTR', start, end, '+',
                                     f"Name={feature_seq.name};ID={feature_id}", feature_id, subject.ref_seq)

    link_sequence_to_feature(feature_seq, feature_rec)


def add_gene_level_subfeature(feature, imgt_feature_name, name_prefix, feature_id, parent_id, ref):
    if imgt_feature_name == 'exon_2':
        name = f"IGHVRegion{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
        add_feature_to_ref(name, 'gene', 'V-REGION', feature['ref_seq'], 'CDS', feature['start']+11, feature['end'], '+',
                           f"Name={feature['gene']}_V-REGION;ID={feature_id}", parent_id, ref)

        name = f"IGHVLpart2{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
        add_feature_to_ref(name, 'gene', 'L_PART-2', feature['ref_seq'], 'CDS', feature['start'], feature['start']+11, '+',
                           f"Name={feature['gene']}_L_PART1;ID={feature_id}", parent_id, ref)
    else:
        name = f"{name_prefix}{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
        add_feature_to_ref(name, 'gene', imgt_feature_name, feature['ref_seq'], 'CDS', feature['start'], feature['end'], '+',
                           f"Name={feature['gene']}_{imgt_feature_name};ID={feature_id}", parent_id, ref)



def add_gene_level_features(session, ref, reference_features):
    feature_id = 1
    for gene, features in reference_features[ref.name].items():
        parent_id = feature_id
        for feature_type, feature in features.items():
            if feature_type == 'gene':
                add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                feature_id += 1
                add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                # TODO - make bed files use IMGT names
            elif feature_type == 'nonamer':
                add_gene_level_subfeature(feature, 'V-NONAMER', 'IGHVNona', feature_id, parent_id, ref)
            elif feature_type == 'spacer':
                add_gene_level_subfeature(feature, 'V-SPACER', 'IGHVSpacer', feature_id, parent_id, ref)
            elif feature_type == 'heptamer':
                add_gene_level_subfeature(feature, 'V-HEPTAMER', 'IGHVHepta', feature_id, parent_id, ref)
                # TODO - split bed exon_2 into L_PART1, V-REGION
            elif feature_type == 'exon_2':
                add_gene_level_subfeature(feature, 'exon_2', 'IGHVRegion', feature_id, parent_id, ref)
            elif feature_type == 'gencode_intron':
                add_gene_level_subfeature(feature, 'V-INTRON', 'IGHVIntron', feature_id, parent_id, ref)
            elif feature_type == 'exon_1':
                add_gene_level_subfeature(feature, 'L_PART-1', 'IGHVLP1', feature_id, parent_id, ref)

        feature_id += 1



