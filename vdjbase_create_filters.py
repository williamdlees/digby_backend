# Create filter definitions for the sample API, using the sqlalchemy class definitions

from db.vdjbase_airr_common import read_definition_data
from db.vdjbase_model import *
from db.vdjbase_airr_model import *

filter_file = 'db/vdjbase_api_query_filters.py'

prelude = '''
# API filter definitions, based on the SQLAlchemy class definitions

# This file is created programmatically by db/vdjbase_create_filters.py. DO NOT UPDATE BY HAND. 

from sqlalchemy import func
from db.vdjbase_airr_model import GenoDetection, Sample, Patient, Study, TissuePro, SeqProtocol, DataPro
from db.vdjbase_model import Allele, Gene, AlleleConfidenceReport

'''

required_classes = {
    'sample_info_filters': [Sample, Patient, Study, TissuePro, SeqProtocol, DataPro, GenoDetection],
    'sequence_filters': [Allele, Gene]
}

extras = {
    'sample_info_filters': """
    'allele': {'model': None, 'field': None},

    'haplotypes': {'model': None, 'field': None},
    'genotypes': {'model': None, 'field': None},

    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset', 'no_uniques': True},
""",
    'sequence_filters': """
    'notes': {'model': 'Allele', 'field': func.group_concat(AlleleConfidenceReport.notes, '\\n').label('notes')},
    'notes_count': {'model': 'Allele', 'field': func.count(AlleleConfidenceReport.id).label('notes_count'), 'sort': 'numeric'},

    'sample_id': {'model': None, 'field': None, 'fieldname': 'sample_id'},
    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset', 'no_uniques': True},
"""
}

excludes = {
    'sample_info_filters': ['id'],
    'sequence_filters': ['id', 'locus_order', 'alpha_order', 'pseudo_gene'],
}

renames = {
    'sample_info_filters': [
        {'name': 'name', 'class': 'Sample', 'rename': 'sample_name'},
        {'name': 'aligner_ver', 'class': 'GenoDetection', 'rename': 'aligner_version'},
        {'name': 'geno_tool', 'class': 'GenoDetection', 'rename': 'Genotyper_Tool'},
        {'name': 'geno_ver', 'class': 'GenoDetection', 'rename': 'Genotyper_Version'},
        {'name': 'haplotype_tool', 'class': 'GenoDetection', 'rename': 'Haplotyper_Tool'},
        {'name': 'haplotype_ver', 'class': 'GenoDetection', 'rename': 'Haplotyper_Version'},
        {'name': 'single_assignment', 'class': 'GenoDetection', 'rename': 'Single Assignment'},
    ],
    'sequence_filters': [
        {'name': 'name', 'class': 'Gene', 'rename': 'gene_name'}
    ],
}



def process_classes(fo, class_defs):
    for filter, class_list in required_classes.items():
        fo.write(f"{filter} = {{\n")

        for required_class in class_list:

            for column in required_class.__table__.columns:
                if column.name in excludes[filter]:
                    continue

                # if column.name == 'aligner_ver':
                #     breakpoint()

                renamed = False
                help_text = ''
                example_text = ''

                if len(column.foreign_keys) > 0:
                    continue

                for rename in renames[filter]:
                    if column.name == rename['name'] and required_class.__name__ == rename['class']:
                        try:
                            renamed = True

                            if required_class.__name__ in class_defs and rename['rename'] in class_defs[required_class.__name__]:
                                help_text = class_defs[required_class.__name__][rename['rename']]['description'].replace("'", '"').replace('\n', '')
                                example_text = class_defs[required_class.__name__][rename['rename']]['example'].replace("'", '"').replace('\n', '')

                            if column.type.python_type is int:
                                fo.write(f"    '{rename['rename']}': {{'model': '{required_class.__name__}', 'field': {required_class.__name__}.{column.name}.label('{rename['rename']}'), 'fieldname': '{column.name}', 'sort': 'numeric', 'help': '{help_text}', 'example': '{example_text}'}},\n")
                            else:
                                fo.write(f"    '{rename['rename']}': {{'model': '{required_class.__name__}', 'field': {required_class.__name__}.{column.name}.label('{rename['rename']}'), 'fieldname': '{column.name}', 'help': '{help_text}', 'example': '{example_text}'}},\n")
                        except NotImplementedError:
                            continue

                if not renamed:
                    try:
                        if required_class.__name__ in class_defs and column.name in class_defs[required_class.__name__]:
                            help_text = class_defs[required_class.__name__][column.name]['description'].replace("'", '"').replace('\n', '')
                            example_text = class_defs[required_class.__name__][column.name]['example'].replace("'", '"').replace('\n', '')
                            
                        if column.type.python_type is int:
                            fo.write(f"    '{column.name}': {{'model': '{required_class.__name__}', 'field': {required_class.__name__}.{column.name}, 'sort': 'numeric', 'help': '{help_text}', 'example': '{example_text}'}},\n")
                        else:
                            fo.write(f"    '{column.name}': {{'model': '{required_class.__name__}', 'field': {required_class.__name__}.{column.name}, 'help': '{help_text}', 'example': '{example_text}'}},\n")
                    except NotImplementedError:
                        continue

        if filter in extras:
            fo.write(extras[filter])

        fo.write('}\n\n')


def write_prelude(fo):
    fo.write(prelude)


def main():
    with open(filter_file, 'w') as fo:
        class_defs = read_definition_data()
        write_prelude(fo)
        process_classes(fo, class_defs)


if __name__ == "__main__":
    main()
