# Rebuild genomic data from files in exports/genomic_metadata

import db.imgt_ref
import os
import yaml
from app import app
from db.digger_assembly import process_digger_record
from db.imgt_assembly import process_imgt_assembly
from db.shared import delete_dependencies
from db.novel_alleles import load_novel_alleles, save_novel_alleles
from db.non_imgt_ref import update_non_imgt_ref

def update_genomic_db():
    status = ""

    #db.imgt_ref.update_imgt()
    status = create_genomic()
    status = status.replace('\n', '<br>')
    return status


def create_genomic():
    genomic_dir = os.path.join(app.config['STATIC_PATH'], 'study_data/Genomic')
    if not os.path.isdir(genomic_dir):
        return 'genomic_metadata folder does not exist - aborting'

    novels = load_novel_alleles(os.path.join(genomic_dir, 'vdjbase_novel_alleles.csv'))
    results = []

    for item in os.listdir(genomic_dir):
        if os.path.isdir(os.path.join(genomic_dir, item)) and item[0] != '.':
            s_dir = os.path.join(genomic_dir, item)
            for s_subdir in os.listdir(s_dir):
                s_subdir = os.path.join(s_dir, s_subdir)
                if os.path.isdir(s_subdir) and s_subdir != '.' and s_subdir != '..':
                    for s_item in os.listdir(s_subdir):
                        s_item_ext = os.path.splitext(s_item)[1]
                        if os.path.isfile(os.path.join(s_subdir, s_item)) and s_item[0] != '.' and s_item_ext in ['.yml', '.yaml']:
                            results.append(create_datasets(os.path.join(s_subdir, s_item), novels))

    save_novel_alleles(novels, os.path.join(genomic_dir, 'vdjbase_novel_alleles.csv'))
    return '\n'.join(results)


def create_datasets(dataset_file, novels):
    results = []

    with open(dataset_file, 'r') as fi:
        study_data = yaml.safe_load(fi)

    ref_sets = {}

    if 'Reference_sets' in study_data:
        for species, file in study_data['Reference_sets'].items():
            ref_sets[species] = os.path.join(os.path.dirname(dataset_file), file)
        del study_data['Reference_sets']

    needed_items = ['Type', 'Annotation_format', 'Species', 'Locus', 'Dataset', 'Sample', 'Institute', 'Researcher', 'Reference', 'Contact', 'Date']

    processed_species = []

    for dataset, sample in study_data.items():
        for sample_name, sample_data in sample.items():
            valid_entry = True
            for needed_item in needed_items:
                if needed_item not in sample_data:
                    results.append('%s: %s not specified' % (dataset_file, needed_item))
                    valid_entry = False
            if not valid_entry:
                continue

            if sample_data['Species'] not in processed_species:
                delete_dependencies(sample_data['Species'])
                if sample_data['Species'] in ref_sets:
                    results.append(update_non_imgt_ref(sample_data['Species'], ref_sets[sample_data['Species']]))
                processed_species.append(sample_data['Species'])

            if sample_data['Annotation_format'] == 'IMGT':
                # results.append(process_imgt_assembly(sample_data))
                pass
            elif sample_data['Annotation_format'] == 'VDJbase':
                results.append(process_digger_record(sample_data, os.path.dirname(dataset_file), novels))
            else:
                results.append('%s: Invalid type/format %s specified for %s'% (dataset_file, sample_data['Type'], sample))

    return "\n".join(results)

