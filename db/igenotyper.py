from hashlib import sha256

import os.path
import shutil
import csv
import glob
from zipfile import ZipFile

from receptor_utils.simple_bio_seq import write_csv

from db.genomic_db import Feature, SampleSequence
from db.genomic_db_functions import save_novel_allele, add_feature_to_ref, find_feature_by_name, \
    find_sequence_by_sequence, save_genomic_sequence, link_sequence_to_feature, update_sample_sequence_link, find_allele_by_name


class GeneParsingException(Exception):
    pass


annotation_records = {}


def read_csv(filename):
    with open(filename, 'r') as fi:
        return list(csv.DictReader(fi))


def to_int(row, field):
    if field in row and row[field]:
        try:
            return int(row[field])
        except ValueError:
            try:
                ret = float(row[field])
                if abs(ret - int(ret)) > 0.0001:
                    raise ValueError(f"{field} is not an integer in {row}")
                return int(ret)
            except ValueError:
                print(f"Error: {field} is not an integer in {row}")
            return 0
    else:
        return 0


def to_float(row, field):
    if field in row and row[field]:
        try:
            return float(row[field])
        except ValueError:
            print(f"Error: {field} is not an float in {row}")
            return 0
    else:
        return 0


def get_gene_type(label):
    gene = label.split('*')[0]

    if gene[3] == 'J' or gene[3] == 'V':
        return gene[3]
    elif gene[3] == 'D' and '-' in gene and '_' not in gene or ('C' not in gene and ('TRD' in gene or 'TRB' in gene)):
        return 'D'
    else:
        return 'C'

    

# TODO - Simplify the schema by merging Feature and Sequence. See ../../notes on genomic schema.txt


def process_genomic_record(session, dataset_dir, sample, annotation_file, reference_features, bam_path):
    global annotation_records

    print(f"Importing sample {sample.sample_name}")

    if annotation_file not in annotation_records:
        annotation_records[annotation_file] = read_csv(annotation_file)

    rows = [x for x in annotation_records[annotation_file] if str(x['sample_name']) == str(sample.sample_id) and str(x['subject']) == str(sample.patient.subject_id) and x['project'] == sample.patient.study.study_id]

    if not os.path.isdir(os.path.join(dataset_dir, 'samples')):
        os.mkdir(os.path.join(dataset_dir, 'samples'))
    
    study_name = sample.patient.study.study_name
    study_path = os.path.join(dataset_dir, 'samples', study_name)

    if not os.path.isdir(study_path):
        os.mkdir(study_path)
        with open(os.path.join(dataset_dir, study_path, 'missing.html'), 'w') as fi:
            fi.write('<html><body><h2>Missing File</h2><p>The file you requested cannot be found. This is likely to be because the annotation pipeline produced no results for this sample.</p></body></html>')

    if not rows:
        # quick fix - try subject as sample name
        rows = [x for x in annotation_records[annotation_file] if str(x['sample_name']) == str(sample.patient.subject_id) and str(x['subject']) == str(sample.patient.subject_id) and x['project'] == sample.patient.study.study_id]

        if not rows:
            print(f'ERROR: {annotation_file} contains no data for sample {sample.sample_id} subject {sample.patient.subject_id} project {sample.patient.study.study_id}')
            sample.annotation_path = '/'.join((study_name, 'missing.html'))
            sample.contig_bam_path = '/'.join((study_name, 'missing.html'))
            return

    # Make the samples directory and project subdirectory if they don't exist

    sample_path = os.path.join(study_path, sample.sample_name)

    if not os.path.isdir(sample_path):
        os.mkdir(sample_path)

    # If the annotation file contains records for multiple samples, split into multiple files

    sample_af_name = f"{sample.sample_name}.csv"

    if len(rows) < len(annotation_records[annotation_file]):
        write_csv(os.path.join(sample_path, sample_af_name), rows)
    else:
        shutil.copy(annotation_file, os.path.join(sample_path, sample_af_name))
    
    sample.annotation_path = '/'.join((study_name, sample.sample_name, sample_af_name))

    # Find the bam files for the sample

    if bam_path:
        if sample.annotation_method == 'IGenotyper':
            bam_path = os.path.join(bam_path, sample.sample_id)
            bamfile_path, bam_files = find_sequence_file(bam_path)
        
            if not bamfile_path:
                sample.contig_bam_path = '/'.join((study_name, 'missing.html'))
                print(f"ERROR: {len(bam_files)} bam files found for sample {sample.sample_name} at path {bam_path}")
            else:
                shutil.copy(bamfile_path, os.path.join(sample_path, f"{sample.sample_name}.bam"))

                if os.path.exists(bamfile_path + '.bai'):
                    shutil.copy(bamfile_path + '.bai', os.path.join(sample_path, f"{sample.sample_name}.bam.bai"))
                else:
                    raise GeneParsingException(f"ERROR: No bai file found for {bamfile_path}")

                sample.contig_bam_path = '/'.join((study_name, sample.sample_name, f"{sample.sample_name}.bam"))
        elif sample.annotation_method == 'Digger':
            fasta_files = glob.glob(os.path.join(bam_path, "*.fasta"))
            fasta_file = [x for x in fasta_files if sample.sample_id in x]
            if len(fasta_file) == 0:
                fasta_file = [x for x in fasta_files if sample.patient.subject_id in x]

            if len(fasta_file) != 1:
                sample.contig_bam_path = '/'.join((study_name, 'missing.html'))
                print(f"ERROR: {len(fasta_file)} fasta files found for sample {sample.sample_name} at path {bam_path}")
            else:
                with ZipFile(os.path.join(sample_path, f"{sample.sample_name}.zip"), 'w') as zip:
                    zip.write(fasta_file[0], f"{sample.sample_name}.fasta")
                sample.contig_bam_path = '/'.join((study_name, sample.sample_name, f"{sample.sample_name}.zip"))

    for row in rows:
        if sample.annotation_method == 'Digger':
            row = convert_digger_row(row)

        if row['vdjbase_allele']:
            seq = process_row(row, session, sample)

            if seq:
                create_features(row, session, sample, reference_features, seq)

    session.commit()


def find_sequence_file(bam_path):
    bamfile_path = None
    bam_files = glob.glob(os.path.join(bam_path, "*.bam"))
    if len(bam_files) == 1:
        bamfile_path = bam_files[0]

    elif len(bam_files) == 0:
        sdir = glob.glob(os.path.join(bam_path, "*"))
        if len(sdir) == 1 and os.path.isdir(sdir[0]):
            bam_files = glob.glob(os.path.join(sdir[0], "*.bam"))
            if len(bam_files) == 1:
                bamfile_path = bam_files[0]

    return bamfile_path, bam_files


def process_row(row, session, sample):
    if 'gene' not in row:
        row['gene'] = row['genotyper_gene']

    if sample.ref_seq and sample.ref_seq.name == 'ighc':
        gene_type = 'C'
    else:
        gene_type = get_gene_type(row['gene'])

    if gene_type in ['V', 'D', 'J', 'C']:
        seq = find_allele_by_name(session, row['vdjbase_allele'])

        if not seq:
            if gene_type == 'V':
                gapped_seq = row['V-REGION-GAPPED']
            else:
                gapped_seq = ''

            seq = save_novel_allele(session, row['gene'], row['vdjbase_allele'], row['notes'].replace('\\n', '\r\n'), row[f'{gene_type}-REGION'], gapped_seq, gene_type)

        update_sample_sequence_link(session, int(row['haplotype'].replace('h=', '')), sample, seq)
        return seq

    return None


def create_features(row, session, sample, reference_features, seq):
    sense = '+'     # by + sense we mean 5' to 3'
    if 'sense' in row:
        sense = row['sense']
    
    if sample.ref_seq and sample.ref_seq.name == 'ighc':
        gene_type = 'C'
    else:
        gene_type = get_gene_type(row['gene'])

    feature = find_feature_by_name(session, f'{gene_type}-REGION', seq.name, sample.ref_seq)

    if feature and seq.sequence.replace('.', '') != feature.feature_seq.replace('.', ''):
        print(f'Error: feature {feature.name} sequence does not match that of sequence {seq.name} in sample {sample.sample_name} subject {sample.patient.patient_name}')

    if gene_type == 'V':
        if not feature:
            feature_id = session.query(Feature).count()

            start = end = 0
            if reference_features:
                if 'REGION' in reference_features[sample.ref_seq.name][seq.gene.name]:
                    start = reference_features[sample.ref_seq.name][seq.gene.name]['REGION']['start']
                    end = reference_features[sample.ref_seq.name][seq.gene.name]['REGION']['end']
                else:
                    start = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_2']['start'] + 11
                    end = reference_features[sample.ref_seq.name][seq.gene.name]['EXON_2']['end']

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

            start = end = 0
            if reference_features:
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
            start = end = 0
            if reference_features:
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
            start = end = 0
            if reference_features:
                ref_seq = sample.ref_seq.name
                start = reference_features[ref_seq][seq.gene.name]['GENE']['start']
                end = reference_features[ref_seq][seq.gene.name]['GENE']['end']

            feature = add_feature_to_ref(seq.name, 'allele', 'C-REGION', seq.sequence, row['C-REGION_CIGAR'], 'CDS', start, end, sense,
                                            f"Name={seq.name}_C-REGION;ID={feature_id}", feature_id, sample.ref_seq)
            link_sequence_to_feature(seq, feature)

            for feature_name, rec in reference_features[ref_seq][seq.gene.name].items():
                if feature_name != 'GENE':
                    add_feature(feature_name, feature_name, reference_features, row, seq, session, sample, sense, ref_seq_name=ref_seq)

    if 'Total_Positions' in row:
        sf = session.query(SampleSequence).filter(SampleSequence.sequence_id == seq.id, SampleSequence.sample_id == sample.id).first()
        sf.total_pos = to_int(row, 'Total_Positions')
        sf.av_coverage = to_float(row, 'Average_Coverage')
        sf.mismatched_positions = to_int(row, 'Mismatched_Positions')
        sf.matched_positions = to_int(row, 'Matched_Positions')
        sf.position_mismatches = row['Position_Mismatches'] if 'Position_Mismatches' in row else 0
        sf.position_matches = row['Position_Matches'] if 'Position_Matches' in row else 0
        sf.percent_accuracy = to_float(row, 'Percent_Accuracy')
        sf.positions_10x = to_int(row, 'Positions_With_At_Least_10x_Coverage')
        sf.fully_spanning_reads = to_int(row, 'Fully_Spanning_Reads')
        sf.fully_spanning_matches = to_int(row, 'Fully_Spanning_Reads_100%_Match')

    # fudge for special ighc bed in human ref
    if sample.ref_seq:
        ref_seq = sample.ref_seq.name
        if gene_type == 'C' and ref_seq == 'igh' and 'ighc' in reference_features:
            ref_seq = 'ighc'
        add_feature('gene_sequence', 'GENE', reference_features, row, seq, session, sample, sense, ref_seq_name=ref_seq)


conv_features = {
    'V': [('V-REGION-GAPPED', 'seq_gapped'), ('V-NONAMER', 'v_nonamer'), ('V-HEPTAMER', 'v_heptamer'), ('V-SPACER', 'spacer_3'), ('L-PART1', 'l_part1'), ('L-PART2', 'l_part2')],
    'D': [('D-REGION', 'seq'), ('D-3_NONAMER', 'd_3_nonamer'), ('D-3_HEPTAMER', 'd_3_heptamer'), ('D-3_SPACER', 'spacer_3'), ('D-5_NONAMER', 'd_5_nonamer'), ('D-5_SPACER', 'spacer_5'), ('D-5_HEPTAMER', 'd_5_heptamer')],
    'J': [('J-REGION', 'seq'), ('J-NONAMER', 'j_nonamer'), ('J-HEPTAMER', 'j_heptamer'), ('J-SPACER', 'spacer_5')],
    'C': [('C-REGION', 'seq')]
}


def convert_digger_row(row):
    row['vdjbase_allele'] = row['asc']
    row['genotyper_gene'] = row['asc_gene']

    if row['genotyper_gene']:
        gene_type = get_gene_type(row['genotyper_gene'])

        row['haplotype'] = '0'

        if gene_type in ['V', 'D', 'J', 'C']:
            row[f'{gene_type}-REGION'] = row['seq']
            row[f'{gene_type}-REGION_CIGAR'] = f'{len(row["seq"])}M'

            for feature, field in conv_features[gene_type]:
                row[feature] = row[field]
                row[feature + '_CIGAR'] = f'{len(row[field])}M'

    return row


# Add a record for a particular sequence observed at a feature if it is not present already. Maintain usage linkages
def add_feature(feature, bed_name, reference_features, row, seq, session, sample, strand='+', ref_seq_name=None):
    if feature not in row or not row[feature]:
        return
    
    if sample.ref_seq and not ref_seq_name:            # fudge for human ighc
        ref_seq_name = sample.ref_seq.name

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

        start = end = 0

        if reference_features:
            start = reference_features[ref_seq_name][seq.gene.name][bed_name]['start']
            end = reference_features[ref_seq_name][seq.gene.name][bed_name]['end']

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
                else:
                    add_gene_level_subfeature(locus, feature, feature_type, f'{locus}CRegion', feature_id, parent_id, ref)
            feature_id += 1


# Add features for igh constant genes. Separate function mainly because of gene naming issues
def add_ighc_gene_level_features(session, ref, reference_features):
    feature_id = 1
    locus = 'IGH'
    for gene, features in reference_features[ref.name].items():
        parent_id = feature_id
        for feature_type, feature in features.items():
            if feature_type == 'GENE':
                add_feature_to_ref(feature['gene'], 'gene', 'C-GENE', feature['ref_seq'], '', 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                feature_id += 1
                add_feature_to_ref(feature['gene'], 'gene', 'C-GENE', feature['ref_seq'], '', 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
            else:
                add_gene_level_subfeature(locus, feature, feature_type, f'{locus}CRegion', feature_id, parent_id, ref)

        feature_id += 1


# Helper function for add_gene_level_features
def add_gene_level_subfeature(locus, feature, imgt_feature_name, name_prefix, feature_id, parent_id, ref):
    name = f"{name_prefix}{feature['gene'][3:]}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
    add_feature_to_ref(name, 'gene', imgt_feature_name, feature['ref_seq'], '', 'CDS', feature['start'], feature['end'], '+',
                        f"Name={feature['gene']}_{imgt_feature_name};ID={feature_id}", parent_id, ref)



