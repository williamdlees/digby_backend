# FUnctions to save genomic elements to the VDJbase database

from db.feature_db import Species, RefSeq, Feature, Sample, SampleSequence, Sequence, Study
from app import db
from Bio.Seq import Seq
from sqlalchemy import and_


def save_genomic_dataset_details(locus, name, ref_file, species, sequence):
    sp = db.session.query(Species).filter_by(name=species).one_or_none()
    if not sp:
        sp = Species(name=species)
        db.session.add(sp)

    ref_seq = RefSeq(name=name, locus=locus, species=sp, sequence=sequence, length=len(sequence))
    db.session.add(ref_seq)


def save_genomic_study(name, institute, researcher, reference, contact, accession_id, accession_reference):
    study = Study(name=name, institute=institute, researcher=researcher, reference=reference, contact=contact, accession_id=accession_id, accession_reference=accession_reference)
    db.session.add(study)
    db.session.commit()
    return study


# Find an allele of this gene that exactly matches the specified sequence
def find_existing_allele(gene_name, gene_sequence):
    gene_sequence = gene_sequence.lower()
    seq = db.session.query(Sequence).filter(and_(Sequence.name.like('%s*i%%' % gene_name), Sequence.sequence == gene_sequence)).one_or_none()

    if seq is None:
        gene_sequence = str(Seq(gene_sequence).reverse_complement()).lower()
    seq = db.session.query(Sequence).filter(and_(Sequence.name.like('%s*i%%' % gene_name), Sequence.sequence == gene_sequence)).one_or_none()

    return seq


# Find all alleles of the specified gene
def find_all_alleles(gene_name):
    sequences = db.session.query(Sequence).filter(Sequence.name.like('%s*i%%' % gene_name)).all()
    return sequences


def save_genomic_sequence(name, imgt_name, novel, deleted, sequence, gapped_sequence, species):
    sequence = Sequence(name=name, imgt_name=imgt_name, type=find_allele_type(name), novel=novel, deleted=deleted, sequence=sequence, gapped_sequence=gapped_sequence, species=species)
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
