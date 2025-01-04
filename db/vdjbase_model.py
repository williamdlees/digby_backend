# coding: utf-8
from sqlalchemy import Boolean, Column, DECIMAL, ForeignKey, Index, Integer, String, Text, func
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from db.db_propertymixin import Details_Mixin

Base = declarative_base()
metadata = Base.metadata


class Details(Details_Mixin, Base):
    __tablename__ = "details"


class Gene(Base):
    __tablename__ = 'gene'

    id = Column(Integer, primary_key=True)
    name = Column(String(250), nullable=False)
    type = Column(String(20), nullable=False)
    family = Column(String(20), nullable=False)
    species = Column(String(20), nullable=False)
    igsnper_plot_path = Column(String(250))
    locus_order = Column(Integer)
    alpha_order = Column(Integer)
    pseudo_gene = Column(Boolean)
    alleles = relationship("Allele")


class HaplotypesFile(Base):
    __tablename__ = 'haplotypes_file'

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


class Allele(Base):
    __tablename__ = 'allele'

    id = Column(Integer, primary_key=True)
    name = Column(String(30), nullable=False)
    pipeline_name = Column(String(30), nullable=False)    # the name assigned by the pipeline, e.g. *bp01
    seq = Column(Text)
    seq_len = Column(String(50), nullable=False)
    similar = Column(String(250), nullable=False)       # list of other names which have the same sequence as this one
    appears = Column(Integer, nullable=False)
    gene_id = Column(ForeignKey('gene.id'), nullable=False, index=True)
    is_single_allele = Column(Boolean, nullable=False)      # False if a short sequence is ambiguous, ie could match more than one reference allele
    low_confidence = Column(Boolean, nullable=False)
    novel = Column(Boolean, nullable=False)
    max_kdiff = Column(DECIMAL, nullable=False)
    closest_ref_id = Column(ForeignKey('allele.id'), index=True)

    closest_ref = relationship('Allele', remote_side=[id])
    gene = relationship('Gene', back_populates="alleles")
    samples = relationship("Sample", secondary="alleles_sample", back_populates="alleles")


class AlleleConfidenceReport(Base):
    __tablename__ = 'allele_confidence_report'

    id = Column(Integer, primary_key=True)
    category = Column(String(50), nullable=False)
    notes = Column(String(400), nullable=False)
    allele_id = Column(ForeignKey('allele.id'), nullable=False, index=True)

    allele = relationship('Allele', backref="confidence")


class AllelesPattern(Base):
    __tablename__ = 'allele_patterns'

    id = Column(Integer, primary_key=True)
    allele_in_p_id = Column(ForeignKey('allele.id'), nullable=False, index=True)
    pattern_id = Column(ForeignKey('allele.id'), nullable=False, index=True)

    allele_in_p = relationship('Allele', primaryjoin='AllelesPattern.allele_in_p_id == Allele.id')
    pattern = relationship('Allele', primaryjoin='AllelesPattern.pattern_id == Allele.id')


class SNP(Base):
    __tablename__ = 'snps'

    id = Column(Integer, primary_key=True)
    pos = Column(Integer, nullable=False)
    from_base = Column(String(1), nullable=False)
    to_base = Column(String(1), nullable=False)
    allele_id = Column(ForeignKey('allele.id'), nullable=False, index=True)

    allele = relationship('Allele')


# The association proxy extension is used for many-many relationshups: see
# https://docs.sqlalchemy.org/en/13/orm/basic_relationships.html

class AllelesSample(Base):
    __tablename__ = 'alleles_sample'

    allele_id = Column(ForeignKey('allele.id'), primary_key=True)
    patient_id = Column(ForeignKey('patient.id'), primary_key=True)
    sample_id = Column(ForeignKey('sample.id'), primary_key=True)
    hap = Column(String(20), nullable=False)
    kdiff = Column(DECIMAL, nullable=False)
    freq_by_clone = Column(Integer, nullable=True)
    freq_by_seq = Column(Integer, nullable=True)
    count = Column(Integer, nullable=True)
    total_count = Column(Integer, nullable=True)

    allele = relationship('Allele', backref='sample_associations', cascade="all")
    patient = relationship('Patient')
    sample = relationship('Sample', backref='allele_associations', cascade="all")


class GenesDistribution(Base):
    __tablename__ = 'genes_distribution'

    id = Column(Integer, primary_key=True)
    frequency = Column(DECIMAL, nullable=False)
    gene_id = Column(ForeignKey('gene.id'), nullable=False, index=True)
    patient_id = Column(ForeignKey('patient.id'), nullable=False, index=True)
    sample_id = Column(ForeignKey('sample.id'), nullable=False, index=True)
    count_by_clones = Column(Boolean, nullable=False)

    gene = relationship('Gene')
    patient = relationship('Patient')
    sample = relationship('Sample')


class HaplotypeEvidence(Base):
    __tablename__ = 'haplotype_evidence'

    id = Column(Integer, primary_key=True)
    allele_id = Column(ForeignKey('allele.id'), nullable=False, index=True)
    sample_id = Column(ForeignKey('sample.id'), nullable=False, index=True)
    hap_gene = Column(String(50), nullable=False)
    counts = Column(String(200), nullable=False)
    scores = Column(String(200), nullable=False)

    allele = relationship('Allele')
    sample = relationship('Sample')


class SamplesHaplotype(Base):
    __tablename__ = 'samples_haplotype'
    __table_args__ = (
        Index('samples_haplotype_samples_id_haplotypes_file_id_50b31fc1_uniq', 'samples_id', 'haplotypes_file_id', unique=True),
    )

    id = Column(Integer, primary_key=True)
    samples_id = Column(ForeignKey('sample.id'), nullable=False, index=True)
    haplotypes_file_id = Column(ForeignKey('haplotypes_file.id'), nullable=False, index=True)

    haplotypes_file = relationship('HaplotypesFile')
    samples = relationship('Sample')

