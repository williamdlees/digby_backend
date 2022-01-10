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
    find_feature_by_sequence, find_sequence_by_sequence, save_genomic_sequence, link_sequence_to_feature, update_subject_sequence_link


class GeneParsingException(Exception):
    pass

annotation_records = {}

def process_igenotyper_record(session, species, dataset_dir, subject, annotation_file, reference_features):
    global annotation_records

    print(f"Importing subject {subject.identifier}")

    if not os.path.isfile(os.path.join(dataset_dir, 'samples')):
        shutil.copy(annotation_file, os.path.join(dataset_dir, 'samples'))

    if annotation_file not in annotation_records:
        annotation_records[annotation_file] = read_eval_genes(annotation_file)

    rows = [x for x in annotation_records[annotation_file] if x['sample_name'] == subject.name_in_study]

    sense = '+'     # by + sense we mean 5' to 3'
    feature_id = 1

    for row in rows:
        chain = row['gene'][3]

        if chain == 'V':
            seq_fields = ['allele_sequence', 'heptamer', 'nonamer', 'spacer']
            # convert record to + sense
            for field in seq_fields:
                row[field] = simple.reverse_complement(row[field])

            # process allele

            seq = find_name_or_novel(session, row)

            if seq is not None:
                update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, seq)
                feature = find_feature_by_name(session, 'V-REGION', seq.name, subject.ref_seq)

                if feature and seq.sequence != feature.feature_seq:
                    print(f'Error: feature {feature.name} sequence does not match that of sequence {seq.name} in subject {subject.identifier}')

                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[subject.ref_seq.name][seq.gene]['gene']['start'] + len(row['leader_sequence'])
                    end = reference_features[subject.ref_seq.name][seq.gene]['gene']['end']
                    feature = add_feature_to_ref(seq.name, 'allele', 'V-REGION', seq.sequence, 'CDS', start, end, '+',
                                   f"Name={seq.name}_V-REGION;ID={feature_id}", feature_id, subject.ref_seq)

                link_sequence_to_feature(seq, feature)

                # process heptamer

                hept_seq = find_sequence_by_sequence(session, 'V-HEPTAMER', seq.gene, row['heptamer'])

                if not hept_seq:
                    hept_name = f"{seq.gene}*{sha256(row['heptamer'].encode('utf-8')).hexdigest()[-4:]}"
                    hept_seq = save_genomic_sequence(session, hept_name, seq.gene, 'V-HEPTAMER', True, False, 'F', row['heptamer'], row['heptamer'])

                update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, hept_seq)
                feature = find_feature_by_name(session, 'V-HEPTAMER', hept_seq.name, subject.ref_seq)

                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[subject.ref_seq.name][seq.gene]['heptamer']['start']
                    end = reference_features[subject.ref_seq.name][seq.gene]['heptamer']['end']
                    feature = add_feature_to_ref(hept_seq.name, 'allele', 'V-HEPTAMER', hept_seq.sequence, 'UTR', start, end, '+',
                                   f"Name={hept_seq.name};ID={feature_id}", feature_id, subject.ref_seq)

                link_sequence_to_feature(hept_seq, feature)

                # process nonamer

                nona_seq = find_sequence_by_sequence(session, 'V-NONAMER', seq.gene, row['nonamer'])

                if not nona_seq:
                    nona_name = f"{seq.gene}*{sha256(row['nonamer'].encode('utf-8')).hexdigest()[-4:]}"
                    nona_seq = save_genomic_sequence(session, nona_name, seq.gene, 'V-NONAMER', True, False, 'F', row['nonamer'], row['nonamer'])

                update_subject_sequence_link(session, int(row['haplotype'].replace('h=', '')), subject, nona_seq)
                feature = find_feature_by_name(session, 'V-NONAMER', nona_seq.name, subject.ref_seq)

                if not feature:
                    feature_id = session.query(Feature).count()
                    start = reference_features[subject.ref_seq.name][seq.gene]['nonamer']['start']
                    end = reference_features[subject.ref_seq.name][seq.gene]['nonamer']['end']
                    feature = add_feature_to_ref(nona_seq.name, 'allele', 'V-NONAMER', nona_seq.sequence, 'UTR', start, end, '+',
                                   f"Name={nona_seq.name};ID={feature_id}", feature_id, subject.ref_seq)

                link_sequence_to_feature(nona_seq, feature)

# For v genes, figure the leader and coding sequence from the entire allele sequence provided by IGenotyper

def find_name_or_novel(session, row):
    try:
        # truncate at the heptamer
        seq = row['allele_sequence']

        if len(seq) < 200:
            raise GeneParsingException(f"Suspiciously short sequence for sample {row['sample_name']} gene {row['gene']} hap {row['haplotype']} (length {len(seq)})")

        heptamer_positions = [_.start() for _ in re.finditer(row['heptamer'], seq)]

        if len(heptamer_positions) == 0:
            raise GeneParsingException(f"heptamer not found in sample {row['sample_name']} gene {row['gene']} hap {row['haplotype']}")

        seq = seq[:heptamer_positions[-1]]

        # find the closest reference allele, using the allele name and gene name provided in the row

        if ref_seq := session.query(Sequence).filter(and_(Sequence.name == row['allele'], Sequence.novel == False)).one_or_none():
            #                    print('allele found in refs')
            allele_seq = ref_seq.sequence
            allele_seq_gapped = ref_seq.gapped_sequence
        elif ref_seq := session.query(Sequence).filter(and_(Sequence.name == row['gene']+'*01', Sequence.novel == False)).one_or_none():
            #                    print('allele %s not found in refs. using allele 01 of specified gene' % row['allele'])
            allele_seq = ref_seq.sequence
            allele_seq_gapped = ref_seq.gapped_sequence
        else:
            seq_choices = session.query(Sequence).filter(and_(Sequence.gene == row['gene'], Sequence.novel == False)).all()
            if not seq_choices:
                raise GeneParsingException('No sequence found in reference set for gene %s (sample %s)' % (row['gene'], row['sample_name']))
            allele_seq = seq_choices[0].sequence
            allele_seq_gapped = seq_choices[0].gapped_sequence
        #                    print('allele *01 not found in refs, and no *01. using %s of specified gene' % (seq_choices[0]))

        # find the best alignment to the reference allele sequence - either an exact match or the closest match
        # within a window close to the end of the gene sequence

        #if row['gene'] == 'IGHV3-30-3':
        #    breakpoint()

        if allele_seq in seq:
            seq_pos = list(re.finditer(allele_seq, seq))[-1]
            row['coding_sequence'] = seq[seq_pos.start():]
            row['leader_sequence'] = seq[:seq_pos.start()]
            row['imgt_allele_match'] = 100.0
            exact_match = True
        else:
            diffs = []
            window = 10
            median_pos = max(0, len(seq) - len(allele_seq))
            best_pos, best_diff = 0, 999

            for i in range(median_pos - window, median_pos + window):
                if i < 0:
                    continue
                if i + len(allele_seq) <= len(seq):
                    diff = simple.nt_diff(seq[i:i + len(allele_seq)], allele_seq)
                else:
                    diff = simple.nt_diff(seq[i:], allele_seq[:len(seq[:len(seq[i:]) - len(allele_seq)])])

                diffs.append(diff)
                if diff < best_diff:
                    best_pos, best_diff = i, diff

            exact_match = False
            row['coding_sequence'] = seq[best_pos:]
            row['leader_sequence'] = seq[:best_pos]
            score = len(row['coding_sequence']) - best_diff
            row['imgt_allele_match'] = round(100 * score / len(row['coding_sequence']), 2) if row['coding_sequence'] else 0

        # assign a name based on closest sequence in the reference set of alleles for this gene
        # If we had an exact match, use just that reference allele - this gives us a good name for extensions

        if not len(row['coding_sequence']):
            raise GeneParsingException(f"Cannot determine gene sequence for sample {row['sample_name']} gene {row['gene']} hap {row['haplotype']}")

        # Now we know the sequence - have we seen it before? If not, store as novel

        s = find_gene_allele_by_seq(session, row['gene'], row['coding_sequence'])

        if not s:
            row['ref_extension'] = False
            if exact_match:
                if len(row['coding_sequence']) > len(allele_seq):
                    row['ref_extension'] = True
                gene_refs = {row['allele']: allele_seq_gapped}
            else:
                seq_choices = session.query(Sequence).filter(and_(Sequence.gene == row['gene'], Sequence.novel == False)).all()
                gene_refs = {seq.name: seq.gapped_sequence for seq in seq_choices}

            name, gapped_sequence = name_novel(row['coding_sequence'], gene_refs, True)
            name, _ = rationalise_name(row['coding_sequence'], name)

            if not exact_match or len(row['coding_sequence']) > len(allele_seq):
                save_novel_allele(session, row['gene'], name, 'F', row['coding_sequence'], gapped_sequence)
        else:
            name = s.name

    except GeneParsingException as e:
        print(e)
        return None

    #print(name)
    s = session.query(Sequence).filter(Sequence.name == name).one_or_none()
    return s

# Read contents of eval_genes.txt into an iterable list, add headers as dict


def read_eval_genes(infile):
    r_headers = ['sample_name',
                 'haplotype',
                 'gene',
                 'allele',
                 'repseq_support',
                 'heptamer',
                 'heptamer_length',
                 'spacer',
                 'spacer_length',
                 'nonamer',
                 'nonamer_length',
                 'allele_sequence',
                 ]

    with open(infile, 'r') as fi:
        reader = csv.DictReader(fi, delimiter='\t', fieldnames=r_headers)
        recs = []
        for row in reader:
            recs.append(row)
        return recs

def add_gene_level_features(session, ref, reference_features):
    feature_id = 1
    for gene, features in reference_features[ref.name].items():
        parent_id = feature_id
        for feature_type, feature in features.items():
            if feature_type == 'gene':
                add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], 'gene', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", feature_id, ref)
                feature_id += 1
                add_feature_to_ref(feature['gene'], 'gene', 'V-GENE', feature['ref_seq'], 'mRNA', feature['start'], feature['end'], '+', f"Name={feature['gene']};ID={feature_id}", parent_id, ref)
            elif feature_type == 'nonamer':
                name = f"IGHVNona{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
                add_feature_to_ref(name, 'gene', 'V-NONAMER', feature['ref_seq'], 'CDS', feature['start'], feature['end'], '+',
                                   f"Name={feature['gene']}_{feature_type};ID={feature_id}", parent_id, ref)
            elif feature_type == 'heptamer':
                name = f"IGHVNona{feature['gene'].replace('IGHV', '')}*{sha256(feature['ref_seq'].encode('utf-8')).hexdigest()[-4:]}"
                add_feature_to_ref(name, 'gene', 'V-HEPTAMER', feature['ref_seq'], 'CDS', feature['start'], feature['end'], '+', f"Name={feature['gene']}_{feature_type};ID={feature_id}", parent_id, ref)
            feature_id += 1



