from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, DateTime
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from db.vdjbase_model import Base


class SeqProtocol(Base):
    __tablename__ = 'database_seq_protocols'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    umi = Column(Boolean, nullable=False)
    read_length = Column(String(20))
    forward_pcr_primer_target_location = Column(String(50))
    reverse_pcr_primer_target_location = Column(String(50))
    sequencing_platform = Column(String(40))
    helix = Column(String(8), nullable=False)


class Study(Base):
    __tablename__ = 'database_studies'

    id = Column(Integer, primary_key=True)
    study_title = Column(String(50), nullable=False)
    lab_address = Column(String(50), nullable=False)
    submitted_by = Column(String(50), nullable=False)
    num_subjects = Column(Integer, nullable=False)
    num_samples = Column(Integer, nullable=False)
    pub_ids = Column(String(200))
    study_contact = Column(String(200))
    study_id = Column(String(50))
    accession_reference = Column(String(200))


class TissuePro(Base):
    __tablename__ = 'database_tissue_pro'

    id = Column(Integer, primary_key=True)
    tissue_processing = Column(String(250), nullable=False)
    cell_species_label = Column(String(20), nullable=False)
    tissue_label = Column(String(100), nullable=False)
    cell_subset_label = Column(String(30), nullable=False)
    sub_cell_type = Column(String(30), nullable=False)
    cell_phenotype = Column(String(30), nullable=False)

    @hybrid_property
    def combined_cell_type(self):
        return self.cell_subset_label + '/' + self.sub_cell_type


class Patient(Base):
    __tablename__ = 'database_patients'

    id = Column(Integer, primary_key=True)
    name = Column(String(50), nullable=False)
    sex = Column(String(20), nullable=False)
    ethnicity = Column(String(50), nullable=False)
    disease_diagnosis_label = Column(String(100), nullable=False)
    study_group_description = Column(String(50))
    age = Column(Integer)
    study_id = Column(ForeignKey('database_studies.id'), nullable=False, index=True)
    name_in_paper = Column(String(30))
    ancestry_population = Column(String(50))
    study = relationship('Study')
    igsnper_sample_id = Column(ForeignKey('database_samples.id'), nullable=False, index=True)
    samples = relationship('Sample', back_populates="patient", primaryjoin="Sample.patient_id==Patient.id")


class Sample(Base):
    __tablename__ = 'database_samples'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    chain = Column(String(30), nullable=False)
    row_reads = Column(Integer, nullable=False)
    genotype = Column(String(100))
    genotype_graph = Column(String(200))
    date = Column(DateTime, nullable=False)
    samples_group = Column(Integer, nullable=False)
    geno_detection_id = Column(ForeignKey('database_geno_detection.id'), nullable=False, index=True)
    patient_id = Column(ForeignKey('database_patients.id'), nullable=False, index=True)
    seq_protocol_id = Column(ForeignKey('database_seq_protocols.id'), nullable=False, index=True)
    study_id = Column(ForeignKey('database_studies.id'), nullable=False, index=True)
    tissue_pro_id = Column(ForeignKey('database_tissue_pro.id'), nullable=False, index=True)
    genotype_stats = Column(String(100))
    genotype_report = Column(String(100))
    igsnper_plot_path = Column(String(250))

    geno_detection = relationship('GenoDetection')
    seq_protocol = relationship('SeqProtocol')
    study = relationship('Study')
    tissue_pro = relationship('TissuePro')
    alleles = relationship("AllelesSample", back_populates="sample")
    patient = relationship('Patient', back_populates='samples', foreign_keys=[patient_id])