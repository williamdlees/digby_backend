from sqlalchemy import Boolean, Column, DECIMAL, DateTime, ForeignKey, Index, Integer, String, Table, Text, func, BigInteger, Float
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class Details(Base):
    __tablename__ = 'details'
    id = Column(Integer, primary_key=True)
    species = Column(String(100))
    locus = Column(String(100))
    date = Column(DateTime())


class RefSeq(Base):
    __tablename__ = 'ref_seq'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    sequence = Column(Text(10000000))
    length = Column(Integer)
    chromosome = Column(String(10))
    start = Column(BigInteger)
    end = Column(BigInteger)
    reference = Column(String(500))
    sense = Column(String(1))
    samples = relationship("Sample", backref='ref_seq')


# The association proxy extension is used for many-many relationships: see
# https://docs.sqlalchemy.org/en/13/orm/basic_relationships.html

class SampleSequence(Base):
    __tablename__ = 'sample_sequence'
    total_pos = Column(Integer, nullable=True)
    av_coverage = Column(Float, nullable=True)
    mismatched_positions = Column(Integer, nullable=True)
    matched_positions = Column(Integer, nullable=True)
    position_mismatches = Column(String, nullable=True)
    position_matches = Column(String, nullable=True)
    percent_accuracy = Column(Float, nullable=True)
    percent_accuracy = Column(Float, nullable=True)
    positions_10x = Column(Integer, nullable=True)
    fully_spanning_reads = Column(Integer, nullable=True)
    fully_spanning_matches = Column(Integer, nullable=True)

    sample_id = Column(ForeignKey('sample.id'), primary_key=True)
    sequence_id = Column(ForeignKey('sequence.id'), primary_key=True)
    sample = relationship('Sample', backref='sequence_associations', cascade="all")
    sequence = relationship('Sequence', backref='sample_associations', cascade="all")
    haplo_count = Column(Integer)
    haplotype = Column(String(50))


class SequenceFeature(Base):
    __tablename__ = 'sequence_feature'
    feature_id = Column(Integer, ForeignKey('feature.id'), primary_key=True)
    sequence_id = Column(Integer, ForeignKey('sequence.id'), primary_key=True)
    feature = relationship('Feature', backref='sequence_associations', cascade="all")
    sequence = relationship('Sequence', backref='feature_associations', cascade="all")


class Feature(Base):
    __tablename__ = 'feature'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    feature_level = Column(String(20))
    feature_type = Column(String(30))
    feature_seq = Column(String(500))
    feature_cigar = Column(String(500))
    feature = Column(String(100))
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String(1))
    attribute = Column(String(200))
    parent_id = Column(Integer)
    refseq_id = Column(Integer, ForeignKey('ref_seq.id'))
    refseq = relationship('RefSeq', backref='features')
    sequences = relationship('Sequence', secondary='sequence_feature', back_populates='features')


class Sequence(Base):
    __tablename__ = 'sequence'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    imgt_name = Column(String(100))
    type = Column(String(100))
    novel = Column(Boolean)
    appearances = Column(Integer)
    deleted = Column(Boolean)
    functional = Column(String(15))
    notes = Column(String(200))
    sequence = Column(Text(1000))
    gapped_sequence = Column(Text(1000))
    max_coverage = Column(Integer)
    max_coverage_sample_id = Column(Integer)
    max_coverage_sample_name = Column(String)
    gene_id = Column(Integer, ForeignKey('gene.id'))
    samples = relationship('Sample', secondary='sample_sequence', back_populates='sequences')
    features = relationship('Feature', secondary='sequence_feature', back_populates='sequences')



class Assembly(Base):
    __tablename__ = 'assembly'
    id = Column(Integer, primary_key=True)
    identifier = Column(String(100))
    reference = Column(String(500))
    sequence_file = Column(String(500))
    sequence = Column(Text(10000000))
    patient_id = Column(Integer, ForeignKey('patient.id'))
    patient = relationship("Patient", backref='assemblies')
    chromosome = Column(String(10))
    start = Column(BigInteger)
    end = Column(BigInteger)


class Gene(Base):
    __tablename__ = 'gene'
    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    type = Column(String(20), nullable=False)
    family = Column(String(20), nullable=False)
    locus_order = Column(Integer)
    alpha_order = Column(Integer)
    pseudo_gene = Column(Boolean)
    sequences = relationship('Sequence', backref='gene')



