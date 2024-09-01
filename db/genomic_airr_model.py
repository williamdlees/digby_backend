
# AIRR-schema related classes in the VDJbase genomic database
# This file is created by vdjbase_create_classes.py. DO NOT UPDATE BY HAND

from sqlalchemy import Boolean, Column, DECIMAL, DateTime, ForeignKey, Index, Integer, String, Table, Text, func, BigInteger, Float
from sqlalchemy.orm import relationship

from db.genomic_db import Base
from db.miairr_mixin import MiAIRR_SampleMixin, MiAIRR_StudyMixin, MiAIRR_PatientMixin, MiAIRR_SeqProtocolMixin, MiAIRR_TissueProMixin, MiAIRR_DataProMixin, MiAIRR_GenoDetectionMixin
from db.db_propertymixin import DB_PropertyMixin


class DB_Properties(DB_PropertyMixin, Base):
    __tablename__ = "db_properties"


class Sample(MiAIRR_SampleMixin, Base):
    __tablename__ = "sample"
    annotation_path = Column(String(100))
    annotation_method = Column(String(100))
    annotation_reference = Column(String(100))
    contig_bam_path = Column(String(100))
    reference_assembly = Column(String(100))

    ref_seq_id = Column(Integer, ForeignKey('ref_seq.id'))
    sequences = relationship('Sequence', secondary='sample_sequence', back_populates='samples')
    patient_id = Column(ForeignKey('patient.id'), nullable=False, index=True)
    seq_protocol_id = Column(ForeignKey('seq_protocol.id'), nullable=False, index=True)
    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    tissue_pro_id = Column(ForeignKey('tissue_pro.id'), nullable=False, index=True)    
    data_pro_id = Column(ForeignKey('data_pro.id'), nullable=False, index=True)    
    seq_protocol = relationship('SeqProtocol')
    study = relationship('Study')
    tissue_pro = relationship('TissuePro')
    data_pro = relationship('DataPro')
    patient = relationship('Patient')


class Study(MiAIRR_StudyMixin, Base):
    __tablename__ = "study"


class Patient(MiAIRR_PatientMixin, Base):
    __tablename__ = "patient"
    mother_in_study = Column(String(100))
    father_in_study = Column(String(100))

    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    study = relationship('Study')
    samples = relationship('Sample', back_populates="patient", primaryjoin="Sample.patient_id==Patient.id")


class TissuePro(MiAIRR_TissueProMixin, Base):
    __tablename__ = "tissue_pro"


class SeqProtocol(MiAIRR_SeqProtocolMixin, Base):
    __tablename__ = "seq_protocol"


class DataPro(MiAIRR_DataProMixin, Base):
    __tablename__ = "data_pro"


class GenoDetection(MiAIRR_GenoDetectionMixin, Base):
    __tablename__ = "geno_detection"

