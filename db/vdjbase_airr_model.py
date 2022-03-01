
# MiAIRR class definitions, based on the MiAIRR schema definition with extensions as noted in vdjbase_airr_schema_defs.xls

# This file is created programmatically by db/vdjbase_create_airr_classes.py. DO NOT UPDATE BY HAND. 


from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, DateTime, DECIMAL
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from db.vdjbase_model import Base



class Study(Base):
    __tablename__ = "study"

    id = Column(Integer, primary_key=True)
    repertoire_id = Column(String(100))
    repertoire_name = Column(String(100))
    repertoire_description = Column(String(100))
    study_id = Column(String(100))
    study_title = Column(String(100))
    study_type_id = Column(String(100))
    study_type_label = Column(String(100))
    study_description = Column(String(100))
    inclusion_exclusion_criteria = Column(String(100))
    grants = Column(String(100))
    study_contact = Column(String(100))
    collected_by = Column(String(100))
    lab_name = Column(String(100))
    lab_address = Column(String(100))
    submitted_by = Column(String(100))
    pub_ids = Column(String(100))
    keywords_study = Column(String(100))
    adc_publish_date = Column(String(100))
    adc_update_date = Column(String(100))
    num_subjects = Column(Integer)
    num_samples = Column(Integer)
    accession_reference = Column(String(100))


class Patient(Base):
    __tablename__ = "patient"

    id = Column(Integer, primary_key=True)
    subject_id = Column(String(100))
    synthetic = Column(Boolean)
    species_id = Column(String(100))
    species_label = Column(String(100))
    organism_id = Column(String(100))
    organism_label = Column(String(100))
    sex = Column(String(100))
    age_min = Column(DECIMAL)
    age_max = Column(DECIMAL)
    age_unit_id = Column(String(100))
    age_unit_label = Column(String(100))
    age_event = Column(String(100))
    age = Column(String(100))
    ancestry_population = Column(String(100))
    ethnicity = Column(String(100))
    race = Column(String(100))
    strain_name = Column(String(100))
    linked_subjects = Column(String(100))
    link_type = Column(String(100))
    study_group_description = Column(String(100))
    disease_diagnosis_id = Column(String(100))
    disease_diagnosis_label = Column(String(100))
    disease_length = Column(String(100))
    disease_stage = Column(String(100))
    prior_therapies = Column(String(100))
    immunogen = Column(String(100))
    intervention = Column(String(100))
    medical_history = Column(String(100))
    receptor_genotype_set_id = Column(String(100))
    mhc_genotype_set_id = Column(String(100))
    mhc_genotype_id = Column(String(100))
    genotype_class = Column(String(100))
    gene_symbol = Column(String(100))
    germline_set_ref = Column(String(100))
    genotype_process = Column(String(100))
    name = Column(String(100))
    igsnper_sample_id = Column(ForeignKey('samples.id'), nullable=False, index=True)

    study_id = Column(ForeignKey('seq_protocol.id'), nullable=False, index=True)
    samples = relationship('Sample', back_populates="patient", primaryjoin="Sample.patient_id==Patient.id")


class Sample(Base):
    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    sample_processing_id = Column(String(100))
    sample_id = Column(String(100))
    sample_type = Column(String(100))
    tissue_id = Column(String(100))
    tissue_label = Column(String(100))
    anatomic_site = Column(String(100))
    disease_state_sample = Column(String(100))
    collection_time_point_relative = Column(DECIMAL)
    collection_time_point_relative_unit_id = Column(String(100))
    collection_time_point_relative_unit_label = Column(String(100))
    collection_time_point_reference = Column(String(100))
    biomaterial_provider = Column(String(100))
    name = Column(String(100))
    reads = Column(Integer)
    genotype = Column(String(100))
    genotype_graph = Column(String(100))
    samples_group = Column(String(100))

    geno_detection_id = Column(ForeignKey('geno_detection.id'), nullable=False, index=True)
    patient_id = Column(ForeignKey('patient.id'), nullable=False, index=True)
    seq_protocol_id = Column(ForeignKey('seq_protocol.id'), nullable=False, index=True)
    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    tissue_pro_id = Column(ForeignKey('tissue_pro.id'), nullable=False, index=True)    
    geno_detection = relationship('GenoDetection')
    seq_protocol = relationship('SeqProtocol')
    study = relationship('Study')
    tissue_pro = relationship('TissuePro')
    alleles = relationship("AllelesSample", back_populates="sample")
    patient = relationship('Patient', back_populates='samples', foreign_keys=[patient_id])


class TissuePro(Base):
    __tablename__ = "tissue_pro"

    id = Column(Integer, primary_key=True)
    tissue_processing = Column(String(100))
    cell_subset_id = Column(String(100))
    cell_subset_label = Column(String(100))
    cell_phenotype = Column(String(100))
    cell_species_id = Column(String(100))
    cell_species_label = Column(String(100))
    single_cell = Column(Boolean)
    cell_number = Column(Integer)
    cells_per_reaction = Column(Integer)
    cell_storage = Column(Boolean)
    cell_quality = Column(String(100))
    cell_isolation = Column(String(100))
    cell_processing_protocol = Column(String(100))
    sub_cell_type = Column(String(100))


class SeqProtocol(Base):
    __tablename__ = "seq_protocol"

    id = Column(Integer, primary_key=True)
    template_class = Column(String(100))
    template_quality = Column(String(100))
    template_amount = Column(DECIMAL)
    template_amount_unit_id = Column(String(100))
    template_amount_unit_label = Column(String(100))
    library_generation_method = Column(String(100))
    library_generation_protocol = Column(String(100))
    library_generation_kit_version = Column(String(100))
    pcr_target_locus = Column(String(100))
    forward_pcr_primer_target_location = Column(String(100))
    reverse_pcr_primer_target_location = Column(String(100))
    complete_sequences = Column(String(100))
    physical_linkage = Column(String(100))
    sequencing_run_id = Column(String(100))
    total_reads_passing_qc_filter = Column(Integer)
    sequencing_platform = Column(String(100))
    sequencing_facility = Column(String(100))
    sequencing_run_date = Column(String(100))
    sequencing_kit = Column(String(100))
    file_type = Column(String(100))
    filename = Column(String(100))
    read_direction = Column(String(100))
    read_length = Column(Integer)
    paired_filename = Column(String(100))
    paired_read_direction = Column(String(100))
    paired_read_length = Column(Integer)
    UMI = Column(Boolean)


class DataPro(Base):
    __tablename__ = "data_pro"

    id = Column(Integer, primary_key=True)
    data_processing_id = Column(String(100))
    primary_annotation = Column(Boolean)
    software_versions = Column(String(100))
    paired_reads_assembly = Column(String(100))
    quality_thresholds = Column(String(100))
    primer_match_cutoffs = Column(String(100))
    collapsing_method = Column(String(100))
    data_processing_protocols = Column(String(100))
    data_processing_files = Column(String(100))
    germline_database = Column(String(100))
    analysis_provenance_id = Column(String(100))

