from app import db
from db.feature_db import Species, RefSeq, Feature, Sample, SampleSequence, Sequence
from Bio import SeqIO
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join, splitext, isdir
import csv
from sqlalchemy import func, and_
from db.shared import delete_dependencies
import traceback

REF_DIR= 'static/pacbio/human/ref'
SAMPLE_DIR = 'static/pacbio/human/samples'



def update():
#    delete_dependencies('Human')
#    return ""

    sp = db.session.query(Species).filter_by(name='Human').one_or_none()

    if not sp:
        sp = Species(name='Human')
        db.session.add(sp)

    for f in listdir(REF_DIR):
        chrom = splitext(f)[0]
        f = join(REF_DIR, f)
        if isfile(f):
            records = list(SeqIO.parse(f, 'fasta'))
            ref_seq = RefSeq(name=chrom, locus='IGH', species=sp, sequence=records[0].seq.lower(), length=len(records[0].seq))
            db.session.add(ref_seq)
            db.session.commit()

    for f in listdir(SAMPLE_DIR):
        p = join(SAMPLE_DIR, f)
        if isdir(p) and f[0] != '.':
            process_sample(p, f)

    return 'Human samples and reference sequences added'


def process_sample(path, sample_name):
    filename = join(path, 'genes_assigned_to_alleles.txt')
    try:
        with open(filename, 'r', newline='') as fi:
            row_count = 0
            feature_id = 1000

            sample = db.session.query(Sample).filter_by(name=sample_name).one_or_none()

            if not sample:
                sample = Sample(name=sample_name, type='Genomic')
                db.session.add(sample)

            gene_id = db.session.query(func.max(Feature.feature_id)).join(RefSeq).join(Species).filter(Species.name == 'Human').filter(Feature.feature == 'gene').count() + 1

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
                        if allele_name != 'ND' and allele_name != 'Deleted':
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
                                            type = 'V-REGION'
                                        elif 'D' in allele_name:
                                            type = 'D-REGION'
                                        elif 'J' in allele_name:
                                            type = 'J-REGION'
                                        else:
                                            type = 'UNKNOWN'

                                        sequence = Sequence(
                                            name=allele_name,
                                            imgt_name='',
                                            type=type,
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
                                    print('Allele name %s not found in IMGT reference set. Ignoring.' % allele_name)
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
                                    type = 'V-REGION'
                                elif 'D' in allele_name:
                                    type = 'D-REGION'
                                elif 'J' in allele_name:
                                    type = 'J-REGION'
                                else:
                                    type = 'UNKNOWN'

                                gene_VDJ = Feature(
                                    name=row['gene_name'] + '_' + type,
                                    feature='CDS',
                                    start=row['start'],
                                    end=row['end'],
                                    strand=strand,
                                    attribute='Name=%s_%s;ID=%04d' % (row['gene_name'], type, gene_id),
                                    feature_id=gene_id,
                                )
                                ref.features.append(gene_VDJ)
                                gene_id += 1

                            sequence.features.append(gene)
                            sequence.features.append(gene_VDJ)
                            db.session.commit()


            except Exception as e:
                print('Unexpected row format: row %d file %s' % (row_count, filename))
                print(e)
                traceback.print_exc()
    except Exception as e:
        print("can't open %s." % filename)
        print(e)

    db.session.commit()
