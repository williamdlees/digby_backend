# AIRR-schema related classes in the VDJbase genomic database

# sample / study / subject

from sqlalchemy import Boolean, Column, DECIMAL, DateTime, ForeignKey, Index, Integer, String, Table, Text, func, BigInteger, Float
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from db.miairr_mixin import MiAIRR_SampleMixin, MiAIRR_StudyMixin, MiAIRR_PatientMixin, MiAIRR_SeqProtocolMixin, MiAIRR_TissueProMixin, MiAIRR_DataProMixin, MiAIRR_GenoDetectionMixin


Base = declarative_base()
metadata = Base.metadata


class Subject(MiAIRR_PatientMixin, Base):
    __tablename__ = 'subject'
#    id = Column(Integer, primary_key=True)     # drop
#    identifier = Column(String(100))            # -> patient_name
#    name_in_study = Column(String(100))         # -> subject_id
    mother_in_study = Column(String(100))       # keep, synthesize
    father_in_study = Column(String(100))       # keep, synthesize
#    age = Column(Integer)                      # drop
#    sex = Column(String(10))                   # drop
#    self_ethnicity = Column(String(100))       # drop
#    grouped_ethnicity = Column(String(100))    # drop
#    population = Column(String(100))           # drop
#    population_abbr = Column(String(100))      # drop
#    super_population = Column(String(100))     # drop

    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    study = relationship('Study')
    samples = relationship('Sample', back_populates="subject", primaryjoin="Sample.subject_id==Subject.id")


class Sample(MiAIRR_SampleMixin, Base):
    __tablename__ = 'sample'
#    identifier = Column(String(100))                # -> patient_name (the VDJbase name)
#    name_in_study = Column(String(100))             # ->sample_name
    annotation_path = Column(String(200))           # keep
#    annotation_method = Column(String(100))         # drop
#    annotation_format = Column(String(100))         # drop
#    annotation_reference = Column(String(100))      # drop - but review: where do we put a ref to Eric's or William's papers?? - in pubids I think
    reference_assembly = Column(String(100))        # keep
#    reference_set_version = Column(String(100))     # ->datapro.germline_database
#    locus_coverage = Column(Float)                  # substitute Eric's coverage metrics
#    sequencing_platform = Column(String(200))       # ->seq_protocol.sequencing_platform
#    assembly_method = Column(String(100))           # drop
#    DNA_source = Column(String(100))                # drop

    ref_seq_id = Column(Integer, ForeignKey('ref_seq.id'))
    sequences = relationship('Sequence', secondary='sample_sequence', back_populates='samples')
    subject_id = Column(ForeignKey('subject.id'), nullable=False, index=True)
    seq_protocol_id = Column(ForeignKey('seq_protocol.id'), nullable=False, index=True)
    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    tissue_pro_id = Column(ForeignKey('tissue_pro.id'), nullable=False, index=True)    
    data_pro_id = Column(ForeignKey('data_pro.id'), nullable=False, index=True)    
    seq_protocol = relationship('SeqProtocol')
    study = relationship('Study')
    tissue_pro = relationship('TissuePro')
    data_pro = relationship('DataPro')
    subject = relationship('Subject')


class Study(MiAIRR_StudyMixin, Base):
    __tablename__ = 'study'
#    study_name = Column(String(50), nullable=False)     # -> study_name
#    study_id = Column(String(50), nullable=False)      # -> study_id
#    title = Column(String(50), nullable=False)         # -> study_title
#    date = Column(DateTime())                          # drop
#    institute = Column(String(500))                    # drop
#    description = Column(String(500))
#    researcher = Column(String(200))                   # drop
#   reference = Column(String(500))                     # drop
#    contact = Column(String(200))                      # drop    
    

class TissuePro(MiAIRR_TissueProMixin, Base):
    __tablename__ = "tissue_pro"


class SeqProtocol(MiAIRR_SeqProtocolMixin, Base):
    __tablename__ = "seq_protocol"


class DataPro(MiAIRR_DataProMixin, Base):
    __tablename__ = "data_pro"



