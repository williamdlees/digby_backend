# Rebuild genomic data from files in exports/genomic_metadata

import db.imgt_ref
import os
import yaml
from app import app
from db.digger_assembly import process_digger_assembly
from db.imgt_assembly import process_imgt_assembly


def update_genomic_db():
    status = ""

    #db.imgt_ref.update_imgt()
    status = create_genomic()
    return status


def create_genomic():
    export_dir = os.path.join(app.config['EXPORT_DIR'], 'genomic_metadata')
    results = []

    if not os.path.isdir(export_dir):
        return 'genomic_metadata folder does not exist - aborting'

    for item in os.listdir(export_dir):
        if os.path.isdir(os.path.join(export_dir, item)) and item[0] != '.':
            s_dir = os.path.join(export_dir, item)
            for s_item in os.listdir(s_dir):
                s_item_ext = os.path.splitext(s_item)[1]
                if os.path.isfile(os.path.join(s_dir, s_item)) and s_item[0] != '.' and s_item_ext in ['.yml', '.yaml']:
                    results.append(create_datasets(os.path.join(s_dir, s_item)))

    return '\n'.join(results)


def create_datasets(dataset_file):
    results = []

    with open(dataset_file, 'r') as fi:
        study_data = yaml.safe_load(fi)

    needed_items = ['Type', 'Species', 'Locus', 'Dataset', 'Sample', 'Institute', 'Researcher', 'Reference', 'Contact', 'Date']

    for dataset, sample in study_data.items():
        for sample_name, sample_data in sample.items():
            valid_entry = True
            for needed_item in needed_items:
                if needed_item not in sample_data:
                    results.append('%s: %s not specified' % (dataset_file, needed_item))
                    valid_entry = False
            if not valid_entry:
                continue

            if sample_data['Type'] == 'IMGT_assembly':
                # results.append(process_imgt_assembly(sample_data))
                pass
            elif sample_data['Type'] == 'VDJbase_assembly':
                results.append(process_digger_assembly(sample_data, os.path.dirname(dataset_file)))
            else:
                results.append('%s: Invalid type %s specified for %s'% (dataset_file, sample_data['Type'], sample))

    return "\n".join(results)

