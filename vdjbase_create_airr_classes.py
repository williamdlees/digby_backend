# Create definitions for the MiAIRR SQLAlchemy classes, using the information from vdjbase_airr_schema_defs.csv

import stringcase

from db.vdjbase_airr_common import read_definition_data

mixin_classes_file = 'db/miairr_mixin.py'
genomic_classes_file = 'db/genomic_airr_model.py'
vdjbase_classes_file = 'db/vdjbase_airr_model.py'

mixin_prelude = '''
# Mixin classes defining MiAIRR fields
# This file is created by vdjbase_create_classes.py. DO NOT UPDATE BY HAND

from sqlalchemy import Column, Integer, String, Boolean, Float


'''

vdjbase_prelude = '''
# MiAIRR class definitions, based on the MiAIRR schema definition with extensions as noted in vdjbase_airr_schema_defs.xls
# This file is created by vdjbase_create_classes.py. DO NOT UPDATE BY HAND


from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, DateTime, Float
from sqlalchemy.orm import relationship

from db.vdjbase_model import Base
from db.miairr_mixin import MiAIRR_SampleMixin, MiAIRR_StudyMixin, MiAIRR_PatientMixin, MiAIRR_SeqProtocolMixin, MiAIRR_TissueProMixin, MiAIRR_DataProMixin, MiAIRR_GenoDetectionMixin

'''

genomic_prelude = '''
# AIRR-schema related classes in the VDJbase genomic database
# This file is created by vdjbase_create_classes.py. DO NOT UPDATE BY HAND

from sqlalchemy import Boolean, Column, DECIMAL, DateTime, ForeignKey, Index, Integer, String, Table, Text, func, BigInteger, Float
from sqlalchemy.orm import relationship

from db.genomic_db import Base
from db.miairr_mixin import MiAIRR_SampleMixin, MiAIRR_StudyMixin, MiAIRR_PatientMixin, MiAIRR_SeqProtocolMixin, MiAIRR_TissueProMixin, MiAIRR_DataProMixin, MiAIRR_GenoDetectionMixin



'''

specials = {
    'igsnper_sample_id': "Column(ForeignKey('sample.id'), nullable=True, index=True)"
}

mixin_key_texts = {}

vdjbase_key_texts = {
    'Patient': """
    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    study = relationship('Study')
    samples = relationship('Sample', back_populates="patient", primaryjoin="Sample.patient_id==Patient.id")
""",

    'Sample': """
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
""",
}

genomic_key_texts = {
    'Patient': """
    study_id = Column(ForeignKey('study.id'), nullable=False, index=True)
    study = relationship('Study')
    samples = relationship('Sample', back_populates="subject", primaryjoin="Sample.subject_id==Subject.id")
""",

    'Sample': """
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
""",
}

def write_prelude(fo, prelude):
    fo.write(prelude)


def write_table(fo, table_name, items, category, key_texts):
    mixin_table_name = f'MiAIRR_{table_name}Mixin'

    if category != 'mixin':
        fo.write(f'''
class {table_name}({mixin_table_name}, Base):
    __tablename__ = "{stringcase.snakecase(table_name)}"
''')
    else:
        fo.write(f'''
class {mixin_table_name}(object):
    id = Column(Integer, primary_key=True)
''')


    for item in items:
        if item['category'] != category:
            continue

        decl = None

        if item['type'] == 'string' or item['type'] == 'genotype' or item['list'] == 'TRUE':
            decl = 'Column(String(100))'
        elif item['type'] == 'number':
            decl = 'Column(Float)'
        elif item['type'] == 'integer':
            decl = 'Column(Integer)'
        elif item['type'] == 'date':
            decl = 'Column(DateTime)'
        elif item['type'] == 'boolean':
            decl = 'Column(Boolean)'

        if not decl:
            if item['type'] == 'haplotype':
                print('Intentionally skipping haplotypes field')
            else:
                print(f'Unrecognised type in {item} - skipped')
            continue

        if item['simple_name'] in specials:
            decl = specials[item['simple_name']]

        fo.write(f"    {item['simple_name']} = {decl}\n")

    if table_name in key_texts:
        fo.write(key_texts[table_name])

    fo.write('\n')




def main():
    defs = read_definition_data()

    with open(mixin_classes_file, 'w', newline='') as fo:
        write_prelude(fo, mixin_prelude)
        for table, items in defs.items():
            write_table(fo, table, items.values(), 'mixin', mixin_key_texts)

    with open(vdjbase_classes_file, 'w', newline='') as fo:
        write_prelude(fo, vdjbase_prelude)
        for table, items in defs.items():
            write_table(fo, table, items.values(), 'airrseq', vdjbase_key_texts)

    with open(genomic_classes_file, 'w', newline='') as fo:
        write_prelude(fo, genomic_prelude)
        for table, items in defs.items():
            write_table(fo, table, items.values(), 'genomic', genomic_key_texts)


if __name__ == "__main__":
    main()
