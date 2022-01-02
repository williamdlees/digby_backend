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


class SubjectSequence(Base):
    __tablename__ = 'subject_sequence'
    id = Column(Integer, primary_key=True, autoincrement=True)
    subject_id = Column(Integer, ForeignKey('subject.id'), nullable=False)
    sequence_id = Column(Integer, ForeignKey('sequence.id'), nullable=False)
    subject = relationship('Subject', backref='sequence_associations', cascade="all")
    sequence = relationship('Sequence', backref='subject_associations', cascade="all")
    haplotype = Column(String(50))


class SequenceFeature(Base):
    __tablename__ = 'sequence_feature'
    id = Column(Integer, primary_key=True, autoincrement=True)
    feature_id = Column(Integer, ForeignKey('feature.id'), nullable=False)
    sequence_id = Column(Integer, ForeignKey('sequence.id'), nullable=False)
    feature = relationship('Feature', backref='sequence_associations', cascade="all")
    sequence = relationship('Sequence', backref='feature_associations', cascade="all")


class Feature(Base):
    __tablename__ = 'feature'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    feature = Column(String(100))
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String(1))
    attribute = Column(String(200))
    parent_id = Column(Integer)
    refseq_id = Column(Integer, ForeignKey('ref_seq.id'))
    refseq = relationship('RefSeq', backref='features')
    sequences = relationship('Sequence', secondary='sequence_feature')


class Sequence(Base):
    __tablename__ = 'sequence'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    gene = Column(String(100))
    imgt_name = Column(String(100))
    type = Column(String(100))
    novel = Column(Boolean)
    deleted = Column(Boolean)
    functional = Column(String(1))
    sequence = Column(Text(10000000))
    gapped_sequence = Column(Text(10000000))
    subject = relationship('Subject', secondary='subject_sequence')
    features = relationship('Feature', secondary='sequence_feature')


class Subject(Base):
    __tablename__ = 'subject'
    id = Column(Integer, primary_key=True)
    identifier = Column(String(100))
    name_in_study = Column(String(100))
    age = Column(Integer)
    sex = Column(String(10))
    annotation_path = Column(String(200))
    annotation_method = Column(String(100))
    annotation_format = Column(String(100))
    annotation_reference = Column(String(100))
    study_id = Column(Integer, ForeignKey('study.id'))
    sequence = relationship('Sequence', secondary='subject_sequence')


class Assembly(Base):
    __tablename__ = 'assembly'
    id = Column(Integer, primary_key=True)
    identifier = Column(String(100))
    reference = Column(String(500))
    sequence_file = Column(String(500))
    sequence = Column(Text(10000000))
    subject_id = Column(Integer, ForeignKey('subject.id'))
    subject = relationship("Subject", backref='assemblies')
    chromosome = Column(String(10))
    start = Column(BigInteger)
    end = Column(BigInteger)


class Study(Base):
    __tablename__ = 'study'
    id = Column(Integer, primary_key=True)
    name = Column(String(50), nullable=False)
    date = Column(DateTime())
    institute = Column(String(500))
    description = Column(String(500))
    researcher = Column(String(200))
    reference = Column(String(500))
    contact = Column(String(200))
    subjects = relationship("Subject", backref='study')





