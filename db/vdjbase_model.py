# coding: utf-8
from sqlalchemy import Boolean, Column, DECIMAL, DateTime, ForeignKey, Index, Integer, String, Table, Text, func
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

Base = declarative_base()
metadata = Base.metadata


class Gene(Base):
    __tablename__ = 'database_gene'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    type = Column(String(20), nullable=False)
    family = Column(String(20), nullable=False)
    species = Column(String(20), nullable=False)
    locus_order = Column(Integer)
    alpha_order = Column(Integer)
    pseudo_gene = Column(Boolean)


class GenoDetection(Base):
    __tablename__ = 'database_geno_detection'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    prepro_tool = Column(String(50), nullable=False)
    aligner_tool = Column(String(50), nullable=False)
    aligner_ver = Column(String(20), nullable=False)
    aligner_reference = Column(String(50))
    geno_tool = Column(String(50), nullable=False)
    geno_ver = Column(String(20), nullable=False)
    haplotype_tool = Column(String(50))
    haplotype_ver = Column(String(20))
    single_assignment = Column(Boolean, nullable=False)
    detection = Column(String(20), nullable=False)


class HaplotypesFile(Base):
    __tablename__ = 'database_haplotypes_files'

    id = Column(Integer, primary_key=True)
    by_gene = Column(String(50), nullable=False)
    allele_col1 = Column(String(50), nullable=False)
    allele_col2 = Column(String(50), nullable=False)
    file = Column(String(100))

    @hybrid_property
    def by_gene_s(self):
        return self.by_gene[4:]

    @by_gene_s.expression
    def by_gene_s(cls):
        return func.substr(cls.by_gene, 4)


class SeqProtocol(Base):
    __tablename__ = 'database_seq_protocols'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    umi = Column(Boolean, nullable=False)
    sequencing_length = Column(String(20))
    primers_3_location = Column(String(50))
    primers_5_location = Column(String(50))
    sequencing_platform = Column(String(40))
    helix = Column(String(8), nullable=False)


class Study(Base):
    __tablename__ = 'database_studies'

    id = Column(Integer, primary_key=True)
    name = Column(String(50), nullable=False)
    institute = Column(String(50), nullable=False)
    researcher = Column(String(50), nullable=False)
    num_subjects = Column(Integer, nullable=False)
    num_samples = Column(Integer, nullable=False)
    reference = Column(String(200))
    contact = Column(String(200))
    accession_id = Column(String(50))
    accession_reference = Column(String(200))


class TissuePro(Base):
    __tablename__ = 'database_tissue_pro'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    species = Column(String(20), nullable=False)
    tissue = Column(String(100), nullable=False)
    cell_type = Column(String(30), nullable=False)
    sub_cell_type = Column(String(30), nullable=False)
    isotype = Column(String(30), nullable=False)

    @hybrid_property
    def combined_cell_type(self):
        return self.cell_type + '/' + self.sub_cell_type


class Allele(Base):
    __tablename__ = 'database_alleles'

    id = Column(Integer, primary_key=True)
    name = Column(String(30), nullable=False)
    seq = Column(Text)
    seq_len = Column(String(50), nullable=False)
    similar = Column(String(250))
    appears = Column(Integer, nullable=False)
    gene_id = Column(ForeignKey('database_gene.id'), nullable=False, index=True)
    is_single_allele = Column(Boolean, nullable=False)
    low_confidence = Column(Boolean, nullable=False)
    novel = Column(Boolean, nullable=False)
    max_kdiff = Column(DECIMAL, nullable=False)
    closest_ref_id = Column(ForeignKey('database_alleles.id'), index=True)

    closest_ref = relationship('Allele', remote_side=[id])
    gene = relationship('Gene')
    samples = relationship("AllelesSample", back_populates="allele")


class Patient(Base):
    __tablename__ = 'database_patients'

    id = Column(Integer, primary_key=True)
    name = Column(String(50), nullable=False)
    sex = Column(String(20), nullable=False)
    ethnic = Column(String(50), nullable=False)
    status = Column(String(100), nullable=False)
    cohort = Column(String(50))
    age = Column(Integer)
    study_id = Column(ForeignKey('database_studies.id'), nullable=False, index=True)
    name_in_paper = Column(String(30))
    country = Column(String(50))

    study = relationship('Study')


class AlleleConfidenceReport(Base):
    __tablename__ = 'database_alleleconfidencereport'

    id = Column(Integer, primary_key=True)
    category = Column(String(50), nullable=False)
    notes = Column(String(400), nullable=False)
    allele_id = Column(ForeignKey('database_alleles.id'), nullable=False, index=True)

    allele = relationship('Allele')


class AllelesPattern(Base):
    __tablename__ = 'database_alleles_patterns'

    id = Column(Integer, primary_key=True)
    allele_in_p_id = Column(ForeignKey('database_alleles.id'), nullable=False, index=True)
    pattern_id = Column(ForeignKey('database_alleles.id'), nullable=False, index=True)

    allele_in_p = relationship('Allele', primaryjoin='AllelesPattern.allele_in_p_id == Allele.id')
    pattern = relationship('Allele', primaryjoin='AllelesPattern.pattern_id == Allele.id')


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

    geno_detection = relationship('GenoDetection')
    patient = relationship('Patient')
    seq_protocol = relationship('SeqProtocol')
    study = relationship('Study')
    tissue_pro = relationship('TissuePro')
    alleles = relationship("AllelesSample", back_populates="sample")


class SNP(Base):
    __tablename__ = 'database_snps'

    id = Column(Integer, primary_key=True)
    pos = Column(Integer, nullable=False)
    from_base = Column(String(1), nullable=False)
    to_base = Column(String(1), nullable=False)
    allele_id = Column(ForeignKey('database_alleles.id'), nullable=False, index=True)

    allele = relationship('Allele')


class AllelesSample(Base):
    __tablename__ = 'database_alleles_samples'

    id = Column(Integer, primary_key=True)
    hap = Column(String(20), nullable=False)
    kdiff = Column(DECIMAL, nullable=False)
    allele_id = Column(ForeignKey('database_alleles.id'), nullable=False, index=True)
    patient_id = Column(ForeignKey('database_patients.id'), nullable=False, index=True)
    sample_id = Column(ForeignKey('database_samples.id'), nullable=False, index=True)

    allele = relationship('Allele', back_populates="samples")
    patient = relationship('Patient')
    sample = relationship('Sample', back_populates="alleles")


class GenesDistribution(Base):
    __tablename__ = 'database_genes_distribution'

    id = Column(Integer, primary_key=True)
    frequency = Column(DECIMAL, nullable=False)
    gene_id = Column(ForeignKey('database_gene.id'), nullable=False, index=True)
    patient_id = Column(ForeignKey('database_patients.id'), nullable=False, index=True)
    sample_id = Column(ForeignKey('database_samples.id'), nullable=False, index=True)
    count_by_clones = Column(Boolean, nullable=False)

    gene = relationship('Gene')
    patient = relationship('Patient')
    sample = relationship('Sample')


class HaplotypeEvidence(Base):
    __tablename__ = 'database_haplotypeevidence'

    id = Column(Integer, primary_key=True)
    allele_id = Column(ForeignKey('database_alleles.id'), nullable=False, index=True)
    sample_id = Column(ForeignKey('database_samples.id'), nullable=False, index=True)
    hap_gene = Column(String(50), nullable=False)
    counts = Column(String(200), nullable=False)

    allele = relationship('Allele')
    sample = relationship('Sample')


class SamplesHaplotype(Base):
    __tablename__ = 'database_samples_haplotype'
    __table_args__ = (
        Index('database_samples_haplotype_samples_id_haplotypes_files_id_50b31fc1_uniq', 'samples_id', 'haplotypes_files_id', unique=True),
    )

    id = Column(Integer, primary_key=True)
    samples_id = Column(ForeignKey('database_samples.id'), nullable=False, index=True)
    haplotypes_files_id = Column(ForeignKey('database_haplotypes_files.id'), nullable=False, index=True)

    haplotypes_files = relationship('HaplotypesFile')
    samples = relationship('Sample')

