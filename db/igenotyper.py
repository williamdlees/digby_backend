from hashlib import sha256

import os.path
import shutil
import csv

from receptor_utils.simple_bio_seq import write_csv

from db.genomic_db import Feature
from db.genomic_db_functions import save_novel_allele, add_feature_to_ref, find_feature_by_name, \
    find_sequence_by_sequence, save_genomic_sequence, link_sequence_to_feature, update_sample_sequence_link, find_allele_by_name


class GeneParsingException(Exception):
    pass


annotation_records = {}


def read_csv(filename):
    with open(filename, 'r') as fi:
        return list(csv.DictReader(fi))


# TODO - Simplify the schema by merging Feature and Sequence. See ../../notes on genomic schema.txt


def process_igenotyper_record(session, dataset_dir, sample, annotation_file, reference_features):
    global annotation_records

    print(f"Importing sample {sample.identifier}")

    if annotation_file not in annotation_records:
        annotation_records[annotation_file] = read_csv(annotation_file)

    rows = [x for x in annotation_records[annotation_file] if str(x['sample_name']) == str(sample.name_in_study) and str(x['subject']) == str(sample.subject.name_in_study) and x['project'] == sample.subject.study.study_id]

    if not rows:
        print(f'ERROR: {annotation_file} contains no data for sample {sample.name_in_study} subject {sample.subject.name_in_study} project {sample.subject.study.study_id}')

    # If the annotation file contains records for multiple samples, split into multiple files

    if len(rows) < len(annotation_records[annotation_file]):
        af_name, af_ext = os.path.splitext(os.path.basename(annotation_file))
        sample_af_name = f"{af_name}_{sample.identifier}{af_ext}"
        write_csv(os.path.join(dataset_dir, 'samples', sample_af_name), rows)
        sample.annotation_path = sample_af_name
    else:
        if not os.path.isfile(os.path.join(dataset_dir, 'samples')):
            shutil.copy(annotation_file, os.path.join(dataset_dir, 'samples'))
        sample.annotation_path = annotation_file

    sense = '+'     # by + sense we mean 5' to 3'
    feature_id = 1

    for row in rows:
        if 'sense' in row:
            sense = row['sense']

        if not row['vdjbase_allele']:
            continue        # probably an incomplete/fragmentary row

        gene_type = row['genotyper_gene'][3]

        if gene_type in ['V', 'D', 'J', 'C']:
            seq = find_allele_by_name(session, row['vdjbase_allele'])

            if not seq:
                if gene_type == 'V':
                    gapped_seq = row['V-REGION-GAPPED']
                else:
                    gapped_seq = ''

                seq = save_novel_allele(session, row['genotyper_gene'], row['vdjbase_allele'], row['notes'].replace('\\n', '\r\n'), row[f'{gene_type}-REGION'], gapped_seq)

            update_sample_sequence_link(session, int(row['haplotype'].replace('h=', '')), sample, seq)
            feature = find_feature_by_name(session, f'{gene_type}-REGION', seq.name, sample.ref_seq)

            if feature and seq.sequence.replace('.', '') != feature.feature_seq.replace('.', ''):
                print(f'Error: feature {feature.name} sequence does not match that of sequence {seq.name} in sample {sample.identifier} subject {sample.subject.identifier}')

            if gene_type == 'V':
                if not feature:
                    feature_id = session.query(Feature).count()

                    if 'REGION' in reference_features[sample.ref_seq.name][seq.gene.name]:
                        start = reference_features[sample.ref_seq.name][seq.gene.name]['REGION']['start']
                        end = reference_features[sample.ref_seq.name][seq.gene.name]['REGION']['end']
                    else:
                        start = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_2']['start'] + 11
                        end = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_2']['end']

                    # if there are gaps in the alignment against the reference assembly, carry these into the feature
                    # TODO - this shouldn't be needed any longer now: gaps are maintained in the CIGAR, not as dots

                    if '.' in row['V-REGION']:
                        feature_seq = row['V-REGION']
                        seq.sequence = row['V-REGION']
                    else:
                        feature_seq = seq.sequence

                    feature = add_feature_to_ref(seq.name, 'allele', 'V-REGION', feature_seq, row['V-REGION_CIGAR'], 'CDS', start, end, sense,
                                                 f"Name={seq.name}_V-REGION;ID={feature_id}", feature_id, sample.ref_seq)

                link_sequence_to_feature(seq, feature)
                add_feature('V-NONAMER', 'NONAMER', reference_features, row, seq, session, sample, sense)
                add_feature('V-SPACER', 'SPACER', reference_features, row, seq, session, sample, sense)
                add_feature('V-HEPTAMER', 'HEPTAMER', reference_features, row, seq, session, sample, sense)
                add_feature('L-PART2', 'L-PART2', reference_features, row, seq, session, sample, sense)
                add_feature('V-INTRON', 'INTRON', reference_features, row, seq, session, sample, sense)
                add_feature('L-PART1', 'EXON_1', reference_features, row, seq, session, sample, sense)
                add_feature('V-UTR', 'UTR', reference_features, row, seq, session, sample, sense)

            elif gene_type == 'D':
                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_1']['start']
                    end = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_1']['end']
                    feature = add_feature_to_ref(seq.name, 'allele', 'D-REGION', seq.sequence, row['D-REGION_CIGAR'], 'CDS', start, end, sense,
                                                 f"Name={seq.name}_D-REGION;ID={feature_id}", feature_id, sample.ref_seq)

                link_sequence_to_feature(seq, feature)
                add_feature('D-3_NONAMER', '3_NONAMER', reference_features, row, seq, session, sample, sense)
                add_feature('D-3_SPACER', '3_SPACER', reference_features, row, seq, session, sample, sense)
                add_feature('D-3_HEPTAMER', '3_HEPTAMER', reference_features, row, seq, session, sample, sense)
                add_feature('D-5_NONAMER', '5_NONAMER', reference_features, row, seq, session, sample, sense)
                add_feature('D-5_SPACER', '5_SPACER', reference_features, row, seq, session, sample, sense)
                add_feature('D-5_HEPTAMER', '5_HEPTAMER', reference_features, row, seq, session, sample, sense)

            elif gene_type == 'J':
                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_1']['start']
                    end = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_1']['end']
                    feature = add_feature_to_ref(seq.name, 'allele', 'J-REGION', seq.sequence, row['J-REGION_CIGAR'], 'CDS', start, end, sense,
                                                 f"Name={seq.name}_J-REGION;ID={feature_id}", feature_id, sample.ref_seq)

                link_sequence_to_feature(seq, feature)
                add_feature('J-NONAMER', 'NONAMER', reference_features, row, seq, session, sample, sense)
                add_feature('J-SPACER', 'SPACER', reference_features, row, seq, session, sample, sense)
                add_feature('J-HEPTAMER', 'HEPTAMER', reference_features, row, seq, session, sample, sense)

            elif gene_type == 'C':
                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[sample.ref_seq.name][seq.gene.name]['GENE']['start']
                    end = reference_features[sample.ref_seq.name][seq.gene.name]['GENE']['end']
                    feature = add_feature_to_ref(seq.name, 'allele', 'C-REGION', seq.sequence, row['C-REGION_CIGAR'], 'CDS', start, end, sense,
                                                 f"Name={seq.name}_C-REGION;ID={feature_id}", feature_id, sample.ref_seq)

                link_sequence_to_feature(seq, feature)

        add_feature('gene_sequence', 'GENE', reference_features, row, seq, session, sample, sense)
    session.commit()


# Add a record for a particular sequence observed at a feature if it is not present already. Maintain usage linkages
def add_feature(feature, bed_name, reference_features, row, seq, session, sample, strand='+'):
    if feature not in row or not row[feature]:
        return

    feature_seq = find_sequence_by_sequence(session, feature, seq.gene.name, row[feature])

    if not feature_seq:
        if feature == 'gene_sequence':
            feature_name = f"{row['vdjbase_allele']}_{sha256(row[feature].encode('utf-8')).hexdigest()[-4:]}"
        else:
            feature_name = f"{seq.gene.name}*{sha256(row[feature].encode('utf-8')).hexdigest()[-4:]}"
        feature_seq = save_genomic_sequence(session, feature_name, seq.gene.name, feature, False, False, '', row[feature], '')

    update_sample_sequence_link(session, int(row['haplotype'].replace('h=', '')), sample, feature_seq)

    feature_rec = find_feature_by_name(session, feature, feature_seq.name, sample.ref_seq)

    if not feature_rec:
        feature_id = session.query(Feature).count()

        #if bed_name == 'gene':
        #    breakpoint()

        start = reference_features[sample.ref_seq.name][seq.gene.name][bed_name]['start']
        end = reference_features[sample.ref_seq.name][seq.gene.name][bed_name]['end']

        feature_rec = add_feature_to_ref(feature_seq.name, 'allele', feature, feature_seq.sequence,  row[feature + '_CIGAR'], 'UTR', start, end, strand,
                                     f"Name={feature_seq.name};ID={feature_id}", feature_id, sample.ref_seq)

    link_sequence_to_feature(feature_seq, feature_rec)


# Add features for each gene defined in the bed files. These are used to annotate the reference track
def add_gene_level_features(session, ref, reference_features):
    feature_id = 1
    for gene, features in reference_features[ref.name].items():
        locus = gene[:3]
        parent_id = feature_id
        for feature_type, feature in features.items():
            if 'V' in feature['gene']:
                if feature_type == 'GENE':
                    add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                elif feature_type == 'NONAMER':
                    add_gene_level_subfeature(locus, feature, 'V-NONAMER', f'{locus}VNona', feature_id, parent_id, ref)
                elif feature_type == 'SPACER':
                    add_gene_level_subfeature(locus, feature, 'V-SPACER', f'{locus}VSpacer', feature_id, parent_id, ref)
                elif feature_type == 'HEPTAMER':
                    add_gene_level_subfeature(locus, feature, 'V-HEPTAMER', f'{locus}VHepta', feature_id, parent_id, ref)
                elif feature_type == 'EXON_2':
                    add_gene_level_subfeature(locus, feature, 'EXON_2', f'{locus}VRegion', feature_id, parent_id, ref)
                elif feature_type == 'INTRON':
                    add_gene_level_subfeature(locus, feature, 'V-INTRON', f'{locus}VIntron', feature_id, parent_id, ref)
                elif feature_type == 'EXON_1':
                    add_gene_level_subfeature(locus, feature, 'L-PART1', f'{locus}VLP1', feature_id, parent_id, ref)
                elif feature_type == 'UTR':
                    add_gene_level_subfeature(locus, feature, 'V-UTR', f'{locus}VUTR', feature_id, parent_id, ref)
                elif feature_type == 'L-PART2':
                    add_gene_level_subfeature(locus, feature, 'L-PART2', f'{locus}VLP2', feature_id, parent_id, ref)
                elif feature_type == 'REGION':
                    add_gene_level_subfeature(locus, feature, 'V-REGION', f'{locus}VRegion', feature_id, parent_id, ref)

            elif 'D' in feature['gene']:
                if feature_type == 'GENE':
                    add_feature_to_ref(feature['gene'], 'gene', 'D-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'D-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                elif feature_type == '3_NONAMER':
                    add_gene_level_subfeature(locus, feature, 'D-3-NONAMER', f'{locus}D3Nona', feature_id, parent_id, ref)
                elif feature_type == '3_SPACER':
                    add_gene_level_subfeature(locus, feature, 'D-3-SPACER', f'{locus}D3Spacer', feature_id, parent_id, ref)
                elif feature_type == '3_HEPTAMER':
                    add_gene_level_subfeature(locus, feature, 'D-3-HEPTAMER', f'{locus}D3Hepta', feature_id, parent_id, ref)
                elif feature_type == '5_NONAMER':
                    add_gene_level_subfeature(locus, feature, 'D-5-NONAMER', f'{locus}D5Nona', feature_id, parent_id, ref)
                elif feature_type == '5_SPACER':
                    add_gene_level_subfeature(locus, feature, 'D-5-SPACER', f'{locus}IGHD5Spacer', feature_id, parent_id, ref)
                elif feature_type == '5_HEPTAMER':
                    add_gene_level_subfeature(locus, feature, 'D-5-HEPTAMER', f'{locus}IGHD5Hepta', feature_id, parent_id, ref)
                elif feature_type == 'EXON_1':
                    add_gene_level_subfeature(locus, feature, 'D-REGION', f'{locus}DRegion', feature_id, parent_id, ref)

            elif 'J' in feature['gene']:
                if feature_type == 'GENE':
                    add_feature_to_ref(feature['gene'], 'gene', 'J-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'J-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
                elif feature_type == 'NONAMER':
                    add_gene_level_subfeature(locus, feature, 'J-NONAMER', f'{locus}JNona', feature_id, parent_id, ref)
                elif feature_type == 'SPACER':
                    add_gene_level_subfeature(locus, feature, 'J-SPACER', f'{locus}JSpacer', feature_id, parent_id, ref)
                elif feature_type == 'HEPTAMER':
                    add_gene_level_subfeature(locus, feature, 'J-HEPTAMER', f'{locus}JHepta', feature_id, parent_id, ref)
                elif feature_type == 'EXON_1':
                    add_gene_level_subfeature(locus, feature, 'J-REGION', f'{locus}JRegion', feature_id, parent_id, ref)

            elif 'C' in feature['gene']:
                if feature_type == 'GENE':
                    add_feature_to_ref(feature['gene'], 'gene', 'C-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                    feature_id += 1
                    add_feature_to_ref(feature['gene'], 'gene', 'C-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
            feature_id += 1


# Helper function for add_gene_level_features
def add_gene_level_subfeature(locus, feature, imgt_feature_name, name_prefix, feature_id, parent_id, ref):
    name = f"{name_prefix}{feature['gene'][3:]}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
    add_feature_to_ref(name, 'gene', imgt_feature_name, feature['ref_seq'], '', 'CDS', feature['start'], feature['end'], '+',
                        f"Name={feature['gene']}_{imgt_feature_name};ID={feature_id}", parent_id, ref)



