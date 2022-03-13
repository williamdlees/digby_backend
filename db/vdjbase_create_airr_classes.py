# Create definitions for the MiAIRR SQLAlchemy classes, using the information from vdjbase_airr_schema_defs.csv

import stringcase

from db.vdjbase_airr_common import read_definition_data

classes_file = 'vdjbase_airr_model.py'

prelude = '''
# MiAIRR class definitions, based on the MiAIRR schema definition with extensions as noted in vdjbase_airr_schema_defs.xls

# This file is created programmatically by db/vdjbase_create_airr_classes.py. DO NOT UPDATE BY HAND. 


from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, DateTime, Float
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from db.vdjbase_model import Base


'''

specials = {
    'igsnper_sample_id': "Column(ForeignKey('sample.id'), nullable=True, index=True)"
}

key_texts = {
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
    alleles = relationship("AllelesSample", back_populates="sample")
    patient = relationship('Patient', back_populates='samples', foreign_keys=[patient_id])
""",


}

def write_prelude(fo):
    fo.write(prelude)


def write_table(fo, table_name, items):
    fo.write(f'''
class {table_name}(Base):
    __tablename__ = "{stringcase.snakecase(table_name)}"

    id = Column(Integer, primary_key=True)
''')


    for item in items:
        decl = None

        if item['type'] == 'string' or item['list'] == 'TRUE':
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

    with open(classes_file, 'w', newline='') as fo:
        write_prelude(fo)
        for table, items in defs.items():
            write_table(fo, table, items.values())


if __name__ == "__main__":
    main()
