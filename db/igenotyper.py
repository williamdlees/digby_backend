from app import db
from db.feature_db import Species, RefSeq, Feature, Sample, SampleSequence, Sequence, Study
from Bio import SeqIO
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join, splitext, isdir
import csv
from sqlalchemy import func, and_
import traceback
from db.save_genomic import save_genomic_dataset_details, save_genomic_study, find_existing_allele, find_all_alleles, save_genomic_sequence, \
    update_sample_sequence_link, add_feature_to_ref, find_allele_type

IGENOTYPER_DIR = 'static/study_data/IGenotyper'


def update():
    update_refs()
    update_studies()
    return 'Human samples and reference sequences added'


def update_refs():
    with open(IGENOTYPER_DIR + '/Ref/ref.txt', 'r') as fo:
        for row in fo:
            (species, locus, name, ref_file) = row.split('\t')
            records = list(SeqIO.parse(IGENOTYPER_DIR + '/Ref/' + ref_file, 'fasta'))
            save_genomic_dataset_details(locus, name, species, sequence=records[0].seq.lower())
            db.session.commit()


def update_studies():
    for study_dir in listdir(IGENOTYPER_DIR):
        if study_dir != 'Ref' and study_dir[0] != '.':
            study_file = join(IGENOTYPER_DIR, study_dir, 'study.tsv')

            sd = {}

            with open(study_file, 'r') as fo:
                for row in fo:
                    k, v = row.split('\t')
                    sd[k] = v.strip()

            study = save_genomic_study(sd['name'], sd['institute'], sd['researcher'], sd['reference'], sd['contact'], sd['accession_id'], sd['accession_reference'])

            for sf in listdir(join(IGENOTYPER_DIR, study_dir)):
                if sf[0] != '.':
                    sample_dir = join(IGENOTYPER_DIR, study_dir, sf)
                    if isdir(sample_dir):
                        process_sample(sample_dir, study)


def process_sample(sample_dir, study):
    try:
        sample_data = {}

        with open(join(sample_dir, 'sample.tsv'), 'r') as fo:
            for row in fo:
                k, v = row.split('\t')
                sample_data[k.strip()] = v.strip()

        species_id = db.session.query(Species.id).filter(Species.name==sample_data['species']).one_or_none()[0]
        ref_id = db.session.query(RefSeq.id).filter(RefSeq.name==sample_data['ref']).one_or_none()[0]

        with open(join(sample_dir, 'genes_assigned_to_alleles.txt'), 'r', newline='') as fi:
            row_count = 0
            feature_id = 1000

            sample = db.session.query(Sample).filter_by(name=sample_data['name'], study=study).one_or_none()

            if not sample:
                report_link = join(sample_dir, 'IGenotyper_report.html')
                if not isfile(report_link):
                    report_link = ''

                sample = Sample(name=sample_data['name'], type='Genomic', date=sample_data['date'], study=study, species_id=species_id, ref_seq_id=ref_id, report_link=report_link)
                db.session.add(sample)

            gene_id = db.session.query(func.max(Feature.feature_id)).join(RefSeq).join(Species).filter(Species.name == 'Human').one_or_none()[0]

            gene_id = 1 if gene_id is None else gene_id+1

            try:
                row_count += 1
                for row in csv.DictReader(fi, delimiter='\t'):
                    #
                    # Hack the chrom name into something we can work with
                    #
                    if row['chrom'] == 'igh':
                        row['chrom'] = 'Human_IGH'

                    strand = '+'

                    row['start'] = int(row['start'])
                    row['end'] = int(row['end'])
                    if row['start'] > row['end']:
                        (row['start'], row['end']) = (row['end'], row['start'])
                        strand = '-'

                    ref = db.session.query(RefSeq).filter(RefSeq.name == row['chrom']).join(Species).filter(Species.name == 'Human').one_or_none()

                    if not ref:
                        print('Reference sequence for locus %s not found.' % row['chrom'])
                        continue

                    if row['end'] > ref.length:
                        print('Row %d Gene %s: end location is past the end of the reference sequence.' % (row_count, row['gene_name']))

                    for h in [1, 2]:
                        allele_name = row['haplotype_%1d_allele' % h]
                        if allele_name != 'ND' and allele_name != 'NotDetermined':
                            if allele_name == 'Novel':
                                gene_sequence = row['haplotype_%1d_novel_sequence' % h]
                                gene_name = row['gene_name']
                                if not len(gene_sequence):
                                    print('Warning: no sequence for novel allele in gene %s' % row['gene name'])

                                sequence = find_existing_allele(gene_name, gene_sequence)

                                if not sequence:
                                    sequences = find_all_alleles(gene_name)

                                    this_ntsequence = str(Seq(gene_sequence).reverse_complement()).lower()
                                    for seq in sequences:
                                        if seq.sequence == this_ntsequence:
                                            sequence = seq

                                    if not sequence:
                                        novel_num = len(sequences) + 1
                                        allele_name = '%s*i%02d' % (row['gene_name'], novel_num)
                                        sequence = save_genomic_sequence(allele_name, '', find_allele_type(allele_name), True, False, this_ntsequence, '', ref.species)
                                update_sample_sequence_link(h, sample, sequence)

                            elif allele_name == 'Deleted':
                                allele_name = '%s_Del' % (row['gene_name'])
                                sequence = db.session.query(Sequence).filter(Sequence.name.like(allele_name)).one_or_none()

                                if not sequence:
                                    sequence = save_genomic_sequence(allele_name, ref, False, True, this_ntsequence,  '', ref.species)

                                update_sample_sequence_link(h, sample, sequence)
                            else:
                                allele_name = '%s*%02d' % (row['gene_name'], int(allele_name))
                                sequence = db.session.query(Sequence).filter(and_(Sequence.imgt_name == allele_name, Sequence.species == ref.species)).one_or_none()
                                if not sequence:
                                    print('Human allele name %s not found in IMGT reference set. Ignoring.' % allele_name)
                                    continue
                                else:
                                    update_sample_sequence_link(h, sample, sequence)

                            gene = db.session.query(Feature).filter(and_(Feature.name == row['gene_name'], Feature.refseq == ref, Feature.feature == 'gene')).one_or_none()
                            allele_type = find_allele_type(row['gene_name'])

                            if not gene:
                                gene = add_feature_to_ref(row['gene_name'], gene, row['start'], row['end'], strand, 'Name=%s;ID=%04d' % (row['gene_name'], gene_id), gene_id, ref)
                                gene_id += 1

                                gene_mRNA = add_feature_to_ref(row['gene_name'], 'mRNA', row['start'], row['end'], strand, 'Name=%s;ID=%04d' % (row['gene_name'], gene_id), gene_id, ref)
                                gene_id += 1

                                gene_VDJ = add_feature_to_ref(row['gene_name'] + '_' + allele_type, 'CDS', row['start'], row['end'], strand, 'Name=%s_%s;ID=%04d' % (row['gene_name'], allele_type, gene_id), gene_id, ref)
                                gene_id += 1
                                sequence.features.append(gene_VDJ)
                            else:
                                gene_VDJ = db.session.query(Feature).filter(and_(Feature.name == row['gene_name'] + '_' + allele_type, Feature.refseq == ref, Feature.feature == 'CDS')).one_or_none()
                                sequence.features.append(gene_VDJ)

                            db.session.commit()


            except Exception as e:
                print('Unexpected row format: row %d file %s' % (row_count, join(sample_dir, 'genes_assigned_to_alleles.txt')))
                print(e)
                traceback.print_exc()
    except Exception as e:
        print("can't open %s." % join(sample_dir, 'genes_assigned_to_alleles.txt'))
        print(e)

    db.session.commit()

