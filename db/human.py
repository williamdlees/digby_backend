from app import db
from db.feature_db import Species, RefSeq, Feature, Sample, SampleFeature, Sequence
from Bio import SeqIO
from os import listdir
from os.path import isfile, join, splitext, isdir
import csv
from sqlalchemy import func, and_
from db.shared import delete_dependencies

REF_DIR= '../pacbio_ig/human/ref'
SAMPLE_DIR = '../pacbio_ig/human/samples'



def update():
    delete_dependencies('Human')
#    return ""

    sp = db.session.query(Species).filter_by(name='Human').one_or_none()

    if not sp:
        sp = Species(name='Human')
        db.session.add(sp)

    for f in listdir(REF_DIR):
        chrom = splitext(f)[0].lower()
        f = join(REF_DIR, f)
        if isfile(f):
            records = list(SeqIO.parse(f, 'fasta'))
            ref_seq = RefSeq(name='Human', locus=chrom, species=sp, sequence=records[0].seq.lower(), length=len(records[0].seq))
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
                    strand = '+'

                    row['start'] = int(row['start'])
                    row['end'] = int(row['end'])
                    if row['start'] > row['end']:
                        (row['start'], row['end']) = (row['end'], row['start'])
                        strand = '-'

                    ref = db.session.query(RefSeq).filter(RefSeq.locus == row['chrom']).join(Species).filter(Species.name == 'Human').one_or_none()

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
                                    novel_num = len(sequences) + 1
                                    allele_name = '%s*i%02d' % (row['gene_name'], novel_num)

                                    if 'V' in allele_name:
                                        type = 'V-GENE'
                                    elif 'D' in allele_name:
                                        type = 'D-GENE'
                                    elif 'J' in allele_name:
                                        type = 'J-GENE'
                                    else:
                                        type = 'UNKNOWN'

                                    sequence = Sequence(
                                        name=allele_name,
                                        imgt_name='',
                                        type=type,
                                        novel=True,
                                        sequence=row['haplotype_%1d_novel_sequence' % h],
                                        gapped_sequence='',
                                        species=ref.species,
                                    )

                                    db.session.add(sequence)
                                    db.session.commit()
                            else:
                                allele_name = '%s*%02d' % (row['gene_name'], int(allele_name))
                                sequence = db.session.query(Sequence).filter(and_(Sequence.imgt_name == allele_name, Sequence.species == ref.species)).one_or_none()
                                if not sequence:
                                    print('Allele name %s not found in IMGT reference set. Ignoring.' % allele_name)
                                    continue

                            gene = db.session.query(Feature).filter(and_(Feature.name == row['gene_name'], Feature.refseq == ref, Feature.feature == 'gene')).one_or_none()

                            if not gene:
                                gene = Feature(
                                    name=row['gene_name'],
                                    feature='gene',
                                    start=row['start'],
                                    end=row['end'],
                                    strand=strand,
                                    attribute='Name=%s;ID=%04d' % (allele_name, gene_id),
                                    feature_id=gene_id,
                                )
                                ref.features.append(gene)
                                gene_id += 1
                                db.session.commit()

                            sf = db.session.query(SampleFeature).filter(and_(SampleFeature.sample == sample, SampleFeature.feature == gene)).one_or_none()
                            if sf:
                                sf.chromosome += ', h%1d' % h
                            else:
                                SampleFeature(sample=sample, feature=gene, chromosome='h%1d' % h)

                            allele = db.session.query(Feature).filter(and_(Feature.name == allele_name, Feature.refseq == ref, Feature.feature == 'mRNA')).one_or_none()

                            if not allele:
                                allele = Feature(
                                    name=allele_name,
                                    feature='CDS',
                                    start=row['start'],
                                    end=row['end'],
                                    strand=strand,
                                    attribute='Name=%s;ID=%04d;Parent=%04d' % (allele_name, feature_id, gene.feature_id),
                                    feature_id=feature_id,
                                    parent_id=gene.feature_id,
                                )
                                feature_id += 1
                                ref.features.append(allele)
                                db.session.commit()

                            sf = db.session.query(SampleFeature).filter(and_(SampleFeature.sample == sample, SampleFeature.feature == allele)).one_or_none()
                            if sf:
                                sf.chromosome += ', h%1d' % h
                            else:
                                SampleFeature(sample=sample, feature=allele, chromosome='h%1d' % h)

                            allele.sequence = sequence

            except Exception as e:
                print('Unexpected row format: row %d file %s' % (row_count, filename))
                print(e)
    except Exception as e:
        print("can't open %s." % filename)
        print(e)

    db.session.commit()
