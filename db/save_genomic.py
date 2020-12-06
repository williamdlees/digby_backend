# FUnctions to save genomic elements to the VDJbase database

from db.feature_db import Species, RefSeq, Feature, Sample, SampleSequence, Sequence, Study, DataSet
from app import db
from Bio.Seq import Seq
from sqlalchemy import and_


feature_type = {
    "5'UTR": 'five_prime_UTR',
    'L-PART1': 'CDS',
    'V-INTRON': 'intron',
    'L-PART2': 'CDS',
    'V-REGION': 'CDS',
    'D-REGION': 'CDS',
    'J-REGION': 'CDS',
    "3'UTR": 'three_prime_UTR',
    'J-HEPTAMER': 'UTR',
    'J-NONAMER': 'UTR',
    'V-HEPTAMER': 'UTR',
    'V-NONAMER': 'UTR',
    "D-3-HEPTAMER": 'UTR',
    "D-£-NONAMER": 'UTR',
    "D-5-HEPTAMER": 'UTR',
    "D-5-NONAMER": 'UTR',
}


def save_genomic_dataset_details(locus, name, species):
    sp = db.session.query(Species).filter_by(name=species).one_or_none()
    if not sp:
        sp = Species(name=species)
        db.session.add(sp)

    data_set = db.session.query(DataSet).filter(DataSet.name == name).filter(DataSet.species == sp).one_or_none()
    if not data_set:
        data_set = DataSet(name=name, locus=locus, species=sp)
        db.session.add(data_set)

    return sp, data_set


def save_genomic_ref_seq(locus, name, sp, ref_sequence, reference, chromosome, start, end):
    ref_seq = RefSeq(name=name, locus=locus, species=sp, sequence=ref_sequence, length=len(ref_sequence), reference=reference, chromosome=chromosome, start=start, end=end)
    db.session.add(ref_seq)
    return ref_seq


def save_genomic_study(name, institute, researcher, reference, contact, description):
    study = Study(name=name, institute=institute, researcher=researcher, reference=reference, contact=contact, description=description)
    db.session.add(study)
    db.session.commit()
    return study


def save_genomic_sample(name, type, date, study, species_id, ref_seq_id, data_set_id, report_link, description, annot_method, annot_ref):
    sample = db.session.query(Sample).filter_by(name=name).one_or_none()

    if not sample:
        sample = Sample(name=name, type=type, date=date, study=study, species_id=species_id, ref_seq_id=ref_seq_id, data_set_id=data_set_id,
                        report_link=report_link, description=description, annot_method=annot_method, annot_ref=annot_ref)
        db.session.add(sample)

    return sample


# Find an allele of this gene that exactly matches the specified sequence
def find_existing_allele(gene_name, gene_sequence):
    gene_sequence = gene_sequence.lower()
    seq = db.session.query(Sequence).filter(and_(Sequence.name.like('%s*i%%' % gene_name), Sequence.sequence == gene_sequence)).one_or_none()

    if seq is None:
        gene_sequence = str(Seq(gene_sequence).reverse_complement()).lower()
    seq = db.session.query(Sequence).filter(and_(Sequence.name.like('%s*i%%' % gene_name), Sequence.sequence == gene_sequence)).one_or_none()

    return seq

def find_allele_by_seq(gene_sequence, species_id):
    gene_sequence = gene_sequence.lower()
    seq = db.session.query(Sequence)\
        .join(Species)\
        .filter(and_(Species.id == species_id, Sequence.sequence == gene_sequence))\
        .all()

    return seq[0] if len(seq) > 0 else None


# Find all alleles of the specified gene
def find_all_alleles(gene_name):
    sequences = db.session.query(Sequence).filter(Sequence.name.like('%s*i%%' % gene_name)).all()
    return sequences


def save_genomic_sequence(name, imgt_name, allele_type, novel, deleted, functional, sequence, gapped_sequence, species):
    sequence = Sequence(name=name, imgt_name=imgt_name, type=allele_type, novel=novel, functional=functional, deleted=deleted, sequence=sequence.lower(), gapped_sequence=gapped_sequence.lower(), species=species)
    db.session.add(sequence)
    return sequence


def update_sample_sequence_link(h, sample, sequence):
    ss = db.session.query(SampleSequence).filter(SampleSequence.sample == sample, SampleSequence.sequence == sequence).one_or_none()
    if ss:
        ss.chromosome = 'h1, h2'
        ss.chromo_count = 2
    else:
        SampleSequence(sample=sample, sequence=sequence, chromosome='h%1d' % h, chromo_count=1)
    db.session.commit()


def add_feature_to_ref(name, feature, start, end,strand, attribute, feature_id, ref):
    gene = Feature(name=name, feature=feature, start=start, end=end, strand=strand, attribute=attribute, feature_id=feature_id)
    ref.features.append(gene)
    return gene


def find_allele_type(allele_name):
    if 'V' in allele_name:
        allele_type = 'V-REGION'
    elif 'D' in allele_name:
        allele_type = 'D-REGION'
    elif 'J' in allele_name:
        allele_type = 'J-REGION'
    else:
        allele_type = 'UNKNOWN'
    return allele_type
