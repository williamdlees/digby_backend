from app import db
from db.feature_db import Species, RefSeq, Feature, Sample, SampleSequence, Sequence, Study
from Bio import SeqIO
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join, splitext, isdir
import csv
from sqlalchemy import func, and_
import traceback


IGENOTYPER_DIR = 'study_data/igenotyper'



def update():
    update_refs()
    update_studies()
    return 'Human samples and reference sequences added'


def update_refs():
    with open(IGENOTYPER_DIR + '/Ref/ref.txt', 'r') as fo:
        for row in fo:
            (species, locus, name, ref_file) = row.split('\t')

            sp = db.session.query(Species).filter_by(name=species).one_or_none()

            if not sp:
                sp = Species(name=species)
                db.session.add(sp)

            records = list(SeqIO.parse(IGENOTYPER_DIR + '/Ref/' + ref_file, 'fasta'))
            ref_seq = RefSeq(name=name, locus=locus, species=sp, sequence=records[0].seq.lower(), length=len(records[0].seq))
            db.session.add(ref_seq)
            db.session.commit()


def update_studies():
    for study_dir in listdir(IGENOTYPER_DIR):
        if study_dir != 'Ref' and study_dir[0] != '.':
            study_file = join(IGENOTYPER_DIR, study_dir, 'study.tsv')

            sd = {}

            with open(study_file, 'r') as fo:
                for row in fo:
                    k, v = row.split('\t')
                    sd[k] = v

            study = Study(name=sd['name'],
                          institute=sd['institute'],
                          researcher=sd['researcher'],
                          reference=sd['reference'],
                          contact=sd['contact'],
                          accession_id=sd['accession_id'],
                          accession_reference=sd['accession_reference'])
            db.session.add(study)
            db.session.commit()

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
                sample = Sample(name=sample_data['name'], type='Genomic', date=sample_data['date'], study=study, species_id=species_id, ref_seq_id=ref_id)
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
                        if allele_name != 'ND' and allele_name != 'NotDetermined' and allele_name != 'Deleted':
                            if allele_name == 'Novel':
                                if not len(row['haplotype_%1d_novel_sequence' % h]):
                                    print('Warning: no sequence for novel allele in gene %s' % row['gene name'])
                                sequence = db.session.query(Sequence).filter(and_(Sequence.name.like('%s*i%%' % row['gene_name']), Sequence.sequence == row['haplotype_%1d_novel_sequence' % h])).one_or_none()

                                if not sequence:
                                    sequences = db.session.query(Sequence).filter(Sequence.name.like('%s*i%%' % row['gene_name'])).all()

                                    this_ntsequence = str(Seq(row['haplotype_%1d_novel_sequence' % h]).reverse_complement()).lower()
                                    for seq in sequences:
                                        if seq.sequence == this_ntsequence:
                                            sequence = seq

                                    if not sequence:
                                        novel_num = len(sequences) + 1
                                        allele_name = '%s*i%02d' % (row['gene_name'], novel_num)

                                        if 'V' in allele_name:
                                            allele_type = 'V-REGION'
                                        elif 'D' in allele_name:
                                            allele_type = 'D-REGION'
                                        elif 'J' in allele_name:
                                            allele_type = 'J-REGION'
                                        else:
                                            allele_type = 'UNKNOWN'

                                        sequence = Sequence(
                                            name=allele_name,
                                            imgt_name='',
                                            type=allele_type,
                                            novel=True,
                                            sequence=this_ntsequence,
                                            gapped_sequence='',
                                            species=ref.species,
                                        )

                                        db.session.add(sequence)
                                        SampleSequence(sample=sample, sequence=sequence, chromosome='h%1d' % h)
                                        db.session.commit()
                                else:
                                    ss = db.session.query(SampleSequence).filter(SampleSequence.sample == sample, SampleSequence.sequence == sequence).one_or_none()
                                    ss.chromosome = 'h1, h2'
                                    db.session.commit()
                            else:
                                allele_name = '%s*%02d' % (row['gene_name'], int(allele_name))
                                sequence = db.session.query(Sequence).filter(and_(Sequence.imgt_name == allele_name, Sequence.species == ref.species)).one_or_none()
                                if not sequence:
                                    print('Human allele name %s not found in IMGT reference set. Ignoring.' % allele_name)
                                    continue
                                else:
                                    ss = db.session.query(SampleSequence).filter(SampleSequence.sample == sample, SampleSequence.sequence == sequence).one_or_none()
                                    if ss:
                                        ss.chromosome = 'h1, h2'
                                    else:
                                        SampleSequence(sample=sample, sequence=sequence, chromosome='h%1d' % h)
                                    db.session.commit()

                            gene = db.session.query(Feature).filter(and_(Feature.name == row['gene_name'], Feature.refseq == ref, Feature.feature == 'gene')).one_or_none()

                            if not gene:
                                gene = Feature(
                                    name=row['gene_name'],
                                    feature='gene',
                                    start=row['start'],
                                    end=row['end'],
                                    strand=strand,
                                    attribute='Name=%s;ID=%04d' % (row['gene_name'], gene_id),
                                    feature_id=gene_id,
                                )
                                ref.features.append(gene)
                                gene_id += 1

                                gene_mRNA = Feature(
                                    name=row['gene_name'],
                                    feature='mRNA',
                                    start=row['start'],
                                    end=row['end'],
                                    strand=strand,
                                    attribute='Name=%s;ID=%04d' % (row['gene_name'], gene_id),
                                    feature_id=gene_id,
                                )
                                ref.features.append(gene_mRNA)
                                gene_id += 1

                                if 'V' in allele_name:
                                    allele_type = 'V-REGION'
                                elif 'D' in allele_name:
                                    allele_type = 'D-REGION'
                                elif 'J' in allele_name:
                                    allele_type = 'J-REGION'
                                else:
                                    allele_type = 'UNKNOWN'

                                gene_VDJ = Feature(
                                    name=row['gene_name'] + '_' + allele_type,
                                    feature='CDS',
                                    start=row['start'],
                                    end=row['end'],
                                    strand=strand,
                                    attribute='Name=%s_%s;ID=%04d' % (row['gene_name'], allele_type, gene_id),
                                    feature_id=gene_id,
                                )
                                ref.features.append(gene_VDJ)
                                gene_id += 1

                            # sequence.features.append(gene)
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
