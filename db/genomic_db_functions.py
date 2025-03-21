# Functions to save genomic elements to the current database

from receptor_utils import simple_bio_seq as simple
from receptor_utils import novel_allele_name
from db.genomic_db import RefSeq, Feature, SampleSequence, Sequence, Details, Assembly, Gene
from db.genomic_airr_model import Sample, Study, Patient, SeqProtocol, TissuePro, DataPro
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
    "D-3-NONAMER": 'UTR',
    "D-5-HEPTAMER": 'UTR',
    "D-5-NONAMER": 'UTR',
}


def save_genomic_dataset_details(session, species, locus, commit_id, branch):
    details = session.query(Details).one_or_none()
    if not details:
        details = Details(
            dbtype='genomic',
            species=species,
            locus=locus,
            created_on=datetime.datetime.now(),
            created_by='digby_backend',
            software_commit_id=commit_id,
            software_branch=branch
            )
        session.add(details)


def save_genomic_ref_seq(session, name, ref_sequence, reference, chromosome, start, end, sense):
    ref_seq = RefSeq(name=name, sequence='', length=len(ref_sequence), reference=reference, chromosome=chromosome, start=start, end=end, sense=sense)
    session.add(ref_seq)
    return ref_seq


def add_feature_to_ref(name, feature_level, feature_type, feature_seq, cigar, feature, start, end, strand, attribute, parent_id, ref):
    gene = Feature(name=name, feature_level=feature_level, feature_type=feature_type, feature_seq=feature_seq, feature_cigar=cigar, feature=feature,
                   start=start, end=end, strand=strand, attribute=attribute, parent_id=parent_id)

    if ref:
        ref.features.append(gene)

    return gene


def save_genomic_sequence(session, name, gene, allele_type, novel, deleted, functional, sequence, gapped_sequence, imgt_name=''):
    gene_id = session.query(Gene.id).filter_by(name=gene).one_or_none()

    if not gene_id:
        print(f"Can't add sequence {name}: gene {gene} not found")
        return None

    gene_id = gene_id[0]

    sequence = Sequence(name=name, gene_id=gene_id, type=allele_type, novel=novel, functional=functional, deleted=deleted, sequence=sequence.upper(),
                        gapped_sequence=gapped_sequence.upper(), imgt_name=imgt_name, appearances=0)
    session.add(sequence)
    return sequence


def save_genomic_patient(identifier, name_in_study, study):
    patient = Patient(identifier=identifier, name_in_study=name_in_study)
    study.patients.append(patient)
    return patient


def save_genomic_sample(session, identifier, patient, name_in_study, reference_assembly, annotation_format):
    sample = Sample(identifier=identifier, name_in_study=name_in_study, annotation_format=annotation_format)
    if reference_assembly:
        ref = session.query(RefSeq).filter(RefSeq.name == reference_assembly).one_or_none()
        if not ref:
            print(f'Error: reference assembly {reference_assembly} not found')
            return None
        #sample.reference_assembly = ref
        ref.samples.append(sample)

    patient.samples.append(sample)
    session.commit()
    return sample


def save_genomic_assembly(identifier, reference, sequence_file, sequence, chromosome, start, end, patient):
    assembly = Assembly(identifier=identifier, reference=reference, sequence_file=sequence_file, sequence=sequence, chromosome=chromosome, start=start, end=end)
    patient.assemblies.append(assembly)
    return assembly


def save_genomic_study(session, name, title, study_id, date, institute, description, researcher, reference, contact):
    study = Study(study_name=name, title=title, study_id=study_id, date=date, institute=institute, description=description, researcher=researcher, reference=reference, contact=contact)
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
    orig_name, gapped_seq, notes = novel_allele_name.name_novel(gene_sequence, ref_set, v_gene)
    name, root_name = rationalise_name(gene_sequence, orig_name)
    gene = session.query(Gene)\
        .join(Sequence, Gene.id == Sequence.gene_id)\
        .filter(Sequence.name == root_name)\
        .one_or_none()

    if not gene:
        print(f'Error: gene {root_name} not found when trying to name novel allele {name}')
        return None

    s = save_novel_allele(session, gene.name, name, notes, gene_sequence, gapped_seq)

    return s


def link_sequence_to_feature(sequence, feature):
    sequence.features.append(feature)
    

def save_novel_allele(session, gene_name, name, notes, sequence, gapped_sequence):
    if not notes:
        functionality = 'Functional'
    elif 'Stop' in notes:
        functionality = 'Pseudogene'
    else:
        functionality = 'ORF'

    gene_id = session.query(Gene.id).filter_by(name=gene_name).one_or_none()

    if not gene_id:
        print(f"Can't add sequence {name}: gene {gene_name} not found")
        return None

    gene_id = gene_id[0]

    s = Sequence(
        name=name,
        gene_id=gene_id,
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
    if name[:2] != 'IG' and name[:2] != 'TR':
        print(f"Unrecognised name format: {name}")
        quit()

    pp = name[:2]
    prefix, name = name.split(pp)
    prefix = prefix + pp
    gene, allele = name.split('*')
    allele_components = allele.split('_')

    root_name = prefix + gene + '*' + allele_components[0]

    if len(allele_components) > 1 and allele_components[1][0] == 'S':
        root_name += '_' + allele_components[1]
        allele_components = allele_components[1:]

    # Avoid infeasibly long names
    if len(allele_components) > 6:
        hash = sha256(gene_sequence.encode('utf-8')).hexdigest()[-5:]
        name = root_name + '_' + hash
    else:
        name = prefix + name
    return name, root_name


# Find all alleles of the specified gene
def find_all_alleles(session, gene_name):
    sequences = session.query(Sequence).filter(Sequence.name.like('%s*i%%' % gene_name)).all()
    return sequences


def update_sample_sequence_link(session, h, sample, sequence):
    ss = session.query(SampleSequence).filter(SampleSequence.sample == sample, SampleSequence.sequence == sequence).one_or_none()

    if not ss:
        sample.sequences.append(sequence)
        sequence.appearances += 1
        session.flush()
        ss = session.query(SampleSequence).filter(SampleSequence.sample == sample, SampleSequence.sequence == sequence).one_or_none()
        ss.haplotype = h
        ss.haplo_count = 1

    else:
        haplo = ss.haplotype.split(',')
        haplo.append('h%1d' % h)
        haplo = ','.join(sorted(haplo))
        ss.haplotype = haplo
        ss.haplo_count += 1

    session.flush()


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


def find_feature_by_sequence(session, feature_type, seq, ref_seq):
    feature = session.query(Feature) \
        .filter(and_(Feature.feature_level == 'allele',
                     Feature.feature_type == feature_type,
                     Feature.feature_seq == seq,
                     Feature.refseq == ref_seq)) \
        .one_or_none()
    return feature


def find_sequence_by_sequence(session, sequence_type, gene, seq):
    seq_object = session.query(Sequence)\
        .join(Gene)\
        .filter(Sequence.gene_id == Gene.id)\
        .filter(Gene.name == gene)\
        .filter(and_(Sequence.type == sequence_type, Sequence.sequence == seq)) \
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


# Once all samples are onboard, calculate the number of subjects that each sequence appears in
def calculate_appearances(session):
    seqs = session.query(Sequence).all()
    for seq in seqs:
        subject_count = session.query(Patient) \
            .join(Sample, Sample.patient_id == Patient.id) \
            .join(SampleSequence, SampleSequence.sample_id == Sample.id) \
            .join(Sequence, Sequence.id == SampleSequence.sequence_id) \
            .filter(Sequence.id == seq.id) \
            .distinct().count()
        seq.appearances = subject_count


# Calculate the sample with the highest coverage for each sequence
def calculate_max_cov_sample(session):
    seqs = session.query(Sequence).all()
    for seq in seqs:
        ss = session.query(SampleSequence).filter(SampleSequence.sequence_id == seq.id).order_by(SampleSequence.fully_spanning_matches.desc()).first()

        if ss is not None:
            seq.max_coverage = ss.fully_spanning_matches
            seq.max_coverage_sample_id = ss.sample_id
            if ss.fully_spanning_matches:
                seq.max_coverage_sample_name = ss.sample.sample_name
