# Functions to save genomic elements to the current database

from receptor_utils import simple_bio_seq as simple
from receptor_utils import novel_allele_name
from db.genomic_db import RefSeq, Feature, Subject, SubjectSequence, Sequence, Study, Details, Assembly
from sqlalchemy import and_
from db.genomic_ref import find_type

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


def save_genomic_dataset_details(session, species, locus):
    details = session.query(Details).one_or_none()
    if not details:
        details = Details(species=species, locus=locus)
        session.add(details)


def save_genomic_ref_seq(session, name, ref_sequence, reference, chromosome, start, end):
    ref_seq = RefSeq(name=name, sequence=ref_sequence, length=len(ref_sequence), reference=reference, chromosome=chromosome, start=start, end=end)
    session.add(ref_seq)
    return ref_seq


def add_feature_to_ref(name, feature, start, end, strand, feature_id, parent_id, ref):
    gene = Feature(name=name, feature=feature, start=start, end=end, strand=strand, feature_id=feature_id, parent_id=parent_id)
    ref.features.append(gene)
    return gene


def save_genomic_sequence(session, name, imgt_name, allele_type, novel, deleted, functional, sequence, gapped_sequence, species):
    sequence = Sequence(name=name, imgt_name=imgt_name, type=allele_type, novel=novel, functional=functional, deleted=deleted, sequence=sequence.lower(), gapped_sequence=gapped_sequence.lower(), species=species)
    session.add(sequence)
    return sequence


def save_genomic_subject(identifier, name_in_study, age, sex, annotation_path, annotation_method, annotation_format, annotation_reference,
                         chromosome, start, end, study):
    subject = Subject(identifier=identifier, name_in_study=name_in_study, age=age, sex=sex, annotation_path=annotation_path, annotation_method=annotation_method,
                      annotation_format=annotation_format, annotation_reference=annotation_reference, chromosome=chromosome, start=start, end=end)
    study.subjects.append(subject)
    return subject


def save_genomic_assembly(identifier, reference, sequence_file, sequence, subject):
    assembly = Assembly(identifier=identifier, reference=reference, sequence_file=sequence_file, sequence=sequence)
    subject.assemblies.append(assembly)
    return assembly


def save_genomic_study(session, name, date, institute, description, researcher, reference, contact):
    study = Study(name=name, date=date, institute=institute, description=description, researcher=researcher, reference=reference, contact=contact)
    session.add(study)
    return study


# Find an allele of this gene that exactly matches the specified sequence
def find_existing_allele(session, gene_name, gene_sequence):
    gene_sequence = gene_sequence.lower()
    seq = session.query(Sequence).filter(and_(Sequence.gene == gene_name, Sequence.sequence == gene_sequence)).one_or_none()

    if seq is None:
        gene_sequence = simple.reverse_complement(gene_sequence)
    seq = session.query(Sequence).filter(and_(Sequence.gene == gene_name, Sequence.sequence == gene_sequence)).one_or_none()

    return seq


def find_allele_by_seq(session, gene_sequence):
    gene_sequence = gene_sequence.upper()
    seq = session.query(Sequence)\
        .filter(Sequence.sequence == gene_sequence)\
        .all()

    return seq[0] if len(seq) > 0 else None


def find_or_assign_allele(session, gene_sequence, v_gene, functional):
    gene_sequence = gene_sequence.lower()
    if s := find_allele_by_seq(session, gene_sequence):
        return s

    ref_set = get_ref_set(session)
    name, gapped_seq = novel_allele_name.name_novel(gene_sequence, ref_set, v_gene)
    root = name.split('_')[0]
    gene = session.query(Sequence.gene).filter(Sequence.name == root).one_or_none()

    s = Sequence(
        name=name,
        gene=gene,
        imgt_name='',
        type=find_type(name),
        sequence=gene_sequence,
        novel=True,
        deleted=False,
        gapped_sequence=gapped_seq,
        functional=functional,
    )
    session.add(s)
    session.commit()

    return s


# Find all alleles of the specified gene
def find_all_alleles(session, gene_name):
    sequences = session.query(Sequence).filter(Sequence.name.like('%s*i%%' % gene_name)).all()
    return sequences


def update_subject_sequence_link(session, h, subject, sequence):
    ss = session.query(SubjectSequence).filter(SubjectSequence.sample == subject, SubjectSequence.sequence == sequence).one_or_none()
    if ss:
        ss.chromosome = 'h1, h2'
        ss.chromo_count = 2
    else:
        SubjectSequence(sample=subject, sequence=sequence, haplotype='h%1d' % h)
    session.commit()


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


def get_ref_set(session):
    ref_set = {}
    refs = session.query(Sequence).filter(Sequence.novel == False).all()

    if refs is not None:
        ref_set = {ref.name: ref.sequence for ref in refs}

    return ref_set
