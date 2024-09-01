
# MiAIRR class definitions, based on the MiAIRR schema definition with extensions as noted in vdjbase_airr_schema_defs.xls
# This file is created by vdjbase_create_classes.py. DO NOT UPDATE BY HAND


from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, DateTime, Float
from sqlalchemy.orm import relationship

from db.vdjbase_model import Base
from db.miairr_mixin import MiAIRR_SampleMixin, MiAIRR_StudyMixin, MiAIRR_PatientMixin, MiAIRR_SeqProtocolMixin, MiAIRR_TissueProMixin, MiAIRR_DataProMixin, MiAIRR_GenoDetectionMixin
from db.db_propertymixin import DB_PropertyMixin


class DB_Properties(DB_PropertyMixin, Base):
    __tablename__ = "db_properties"


class Sample(MiAIRR_SampleMixin, Base):
    __tablename__ = "sample"
    reads = Column(Integer)
    genotype = Column(String(100))
    igsnper_plot_path = Column(String(100))
    sample_group = Column(String(100))
    genotype_stats = Column(String(100))
    genotype_report = Column(String(100))

    geno_detection_id = Column(ForeignKey('geno_detection.id'), nullable=True, index=True)
    patient_id = Column(ForeignKey('patient.id'), nullable=False, index=True)
    seq_protocol_id = Column(ForeignKey('seq_protocol.id'), nullable=False, index=True)
    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    tissue_pro_id = Column(ForeignKey('tissue_pro.id'), nullable=False, index=True)    
    data_pro_id = Column(ForeignKey('data_pro.id'), nullable=False, index=True)    
    geno_detection = relationship('GenoDetection')
    seq_protocol = relationship('SeqProtocol')
    study = relationship('Study')
    tissue_pro = relationship('TissuePro')
    data_pro = relationship('DataPro')
    alleles = relationship("Allele", secondary="alleles_sample", back_populates="samples")
    patient = relationship('Patient', back_populates='samples', foreign_keys=[patient_id])


class Study(MiAIRR_StudyMixin, Base):
    __tablename__ = "study"


class Patient(MiAIRR_PatientMixin, Base):
    __tablename__ = "patient"
    igsnper_sample_id = Column(ForeignKey('sample.id'), nullable=True, index=True)

    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    study = relationship('Study')
    samples = relationship('Sample', back_populates="patient", primaryjoin="Sample.patient_id==Patient.id")


class TissuePro(MiAIRR_TissueProMixin, Base):
    __tablename__ = "tissue_pro"
    sub_cell_type = Column(String(100))


class SeqProtocol(MiAIRR_SeqProtocolMixin, Base):
    __tablename__ = "seq_protocol"


class DataPro(MiAIRR_DataProMixin, Base):
    __tablename__ = "data_pro"


class GenoDetection(MiAIRR_GenoDetectionMixin, Base):
    __tablename__ = "geno_detection"
    prepro_tool = Column(String(100))
    aligner_tool = Column(String(100))
    aligner_ver = Column(String(100))
    aligner_reference_v = Column(String(100))
    aligner_reference_d = Column(String(100))
    aligner_reference_j = Column(String(100))
    geno_tool = Column(String(100))
    geno_ver = Column(String(100))
    haplotype_tool = Column(String(100))
    haplotype_ver = Column(String(100))
    single_assignment = Column(Boolean)
    sample_basis = Column(String(100))

