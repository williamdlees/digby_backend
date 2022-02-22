# Functions to save genomic elements to the current database

from receptor_utils import simple_bio_seq as simple
from receptor_utils import novel_allele_name
from db.genomic_db import RefSeq, Feature, Subject, SubjectSequence, Sequence, Study, Details, Assembly, SequenceFeature, Gene
from sqlalchemy import and_
from db.genomic_ref import find_type
import datetime
from hashlib import sha256

gff_feature_type = {
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
        details = Details(species=species, locus=locus, date=datetime.datetime.now())
        session.add(details)


def save_genomic_ref_seq(session, name, ref_sequence, reference, chromosome, start, end):
    ref_seq = RefSeq(name=name, sequence=ref_sequence, length=len(ref_sequence), reference=reference, chromosome=chromosome, start=start, end=end)
    session.add(ref_seq)
    return ref_seq


def add_feature_to_ref(name, feature_level, feature_type, feature_seq, feature, start, end, strand, attribute, parent_id, ref):
    gene = Feature(name=name, feature_level=feature_level, feature_type=feature_type, feature_seq=feature_seq, feature=feature,
                   start=start, end=end, strand=strand, attribute=attribute, parent_id=parent_id)
    ref.features.append(gene)
    return gene


def save_genomic_sequence(session, name, gene, allele_type, novel, deleted, functional, sequence, gapped_sequence, imgt_name=''):
    sequence = Sequence(name=name, gene=gene, type=allele_type, novel=novel, functional=functional, deleted=deleted, sequence=sequence.upper(),
                        gapped_sequence=gapped_sequence.upper(), imgt_name=imgt_name, appearances=0)
    session.add(sequence)
    return sequence


def save_genomic_subject(session, identifier, name_in_study, age, sex, annotation_path, annotation_method, annotation_format, annotation_reference, reference_assembly, study):
    subject = Subject(identifier=identifier, name_in_study=name_in_study, age=age, sex=sex, annotation_path=annotation_path, annotation_method=annotation_method,
                      annotation_format=annotation_format, annotation_reference=annotation_reference)
    if reference_assembly:
        ref = session.query(RefSeq).filter(RefSeq.name == reference_assembly).one_or_none()
        if not ref:
            print(f'Error: reference assembly {reference_assembly} not found')
            return None
        ref.subjects.append(subject)

    study.subjects.append(subject)
    return subject


def save_genomic_assembly(identifier, reference, sequence_file, sequence, chromosome, start, end, subject):
    assembly = Assembly(identifier=identifier, reference=reference, sequence_file=sequence_file, sequence=sequence, chromosome=chromosome, start=start, end=end)
    subject.assemblies.append(assembly)
    return assembly


def save_genomic_study(session, name, date, institute, description, researcher, reference, contact):
    study = Study(name=name, date=date, institute=institute, description=description, researcher=researcher, reference=reference, contact=contact)
    session.add(study)
    return study


# Find an allele of this gene that exactly matches the specified sequence
def find_existing_allele(session, gene_name, gene_sequence):
    gene_sequence = gene_sequence.upper()
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


def find_allele_by_name(session, allele_name):
    seq = session.query(Sequence)\
        .filter(Sequence.name == allele_name)\
        .one_or_none()

    return seq


def find_gene_allele_by_seq(session, gene, gene_sequence):
    gene_sequence = gene_sequence.upper()
    seq = session.query(Sequence)\
        .filter(Sequence.sequence == gene_sequence)\
        .filter(Sequence.gene == gene)\
        .all()

    return seq[0] if len(seq) > 0 else None


def find_or_assign_allele(session, gene_sequence, v_gene, functional):
    if s := find_allele_by_seq(session, gene_sequence):
        return s

    ref_set = get_ref_set(session)
    name, gapped_seq, notes = novel_allele_name.name_novel(gene_sequence, ref_set, v_gene)
    name, root_name = rationalise_name(gene_sequence, name)
    gene = session.query(Sequence.gene).filter(Sequence.name == root_name).one_or_none()[0]

    s = save_novel_allele(session, gene, name, notes, gene_sequence, gapped_seq)

    return s


def link_sequence_to_feature(sequence, feature):
    sf = SequenceFeature()
    sf.sequence = sequence
    sf.feature = feature


def save_novel_allele(session, gene, name, notes, sequence, gapped_sequence):
    if not notes:
        functionality = 'Functional'
    elif 'Stop' in notes:
        functionality = 'Pseudogene'
    else:
        functionality = 'ORF'

    s = Sequence(
        name=name,
        gene=gene,
        imgt_name='',
        type=find_type(name),
        sequence=sequence,
        novel=True,
        appearances=0,
        deleted=False,
        gapped_sequence=gapped_sequence,
        functional=functionality,
        notes=notes,
    )
    session.add(s)
    return s


def rationalise_name(gene_sequence, name):
    prefix, name = name.split('IG')
    prefix = prefix + 'IG'
    name_components = name.split('_')
    root_name = prefix + name_components[0]
    # Avoid infeasibly long names
    if len(name_components) > 6:
        hash = sha256(gene_sequence.encode('utf-8')).hexdigest()[-5:]
        name = root_name + '_' + hash
    else:
        name = prefix + name
    return name, root_name


# Find all alleles of the specified gene
def find_all_alleles(session, gene_name):
    sequences = session.query(Sequence).filter(Sequence.name.like('%s*i%%' % gene_name)).all()
    return sequences


def update_subject_sequence_link(session, h, subject, sequence):
    ss = session.query(SubjectSequence).filter(SubjectSequence.subject == subject, SubjectSequence.sequence == sequence).one_or_none()
    if ss:
        haplo = ss.haplotype.split(',')
        haplo.append('h%1d' % h)
        haplo = ','.join(sorted(haplo))
        ss.haplotype = haplo
        ss.haplo_count += 1
    else:
        SubjectSequence(subject=subject, sequence=sequence, haplotype='h%1d' % h, haplo_count=1)
        sequence.appearances += 1
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
        ref_set = {ref.study_title: ref.sequence for ref in refs}

    return ref_set


def find_feature_by_sequence(session, feature_type, seq, ref_seq):
    feature = session.query(Feature) \
        .filter(and_(Feature.feature_level == 'allele',
                     Feature.feature_type == feature_type,
                     Feature.feature_seq == seq,
                     Feature.refseq == ref_seq)) \
        .one_or_none()
    return feature


def find_sequence_by_sequence(session, sequence_type, gene, seq):
    seq_object = session.query(Sequence) \
        .filter(and_(Sequence.type == sequence_type, Sequence.sequence == seq, Sequence.gene == gene)) \
        .one_or_none()
    return seq_object


def find_feature_by_name(session, feature_type, name, ref_seq):
    feature = session.query(Feature) \
        .filter(and_(Feature.feature_level == 'allele',
                     Feature.feature_type == feature_type,
                     Feature.name == name,
                     Feature.refseq == ref_seq)) \
        .one_or_none()
    return feature


