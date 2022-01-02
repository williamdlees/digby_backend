# Create database and associated files for a single genomic dataset
from receptor_utils import simple_bio_seq as simple
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool

import os
import yaml

from db.build_gff import build_gff
from db.digger_assembly import process_digger_record
from db.genomic_db import Base
from db.genomic_db_functions import save_genomic_study, save_genomic_subject, save_genomic_assembly, save_genomic_dataset_details
#from db.imgt_assembly import process_imgt_assembly
from db.genomic_ref import update_genomic_ref


Session = sessionmaker()


def create_dataset(dataset_dir):
    dataset_dir = os.path.abspath(dataset_dir)

    if not os.path.isdir(dataset_dir):
        print(f'Directory {dataset_dir} does not exist.')
        return

    yml_files = []
    for entry in os.scandir(dataset_dir):
        if entry.is_file() and os.path.splitext(entry.name)[1] in ['.yml', '.yaml']:
            yml_files.append(entry.name)

    if len(yml_files) == 0:
        print(f'Error: no yml file found in directory {dataset_dir}.')
        return
    elif len(yml_files) > 1:
        print(f'Error: multiple yml files found in directory {dataset_dir}.')
        return

    dataset_file = yml_files[0]
    with open(dataset_file, 'r') as fi:
        study_data = yaml.safe_load(fi)

    dataset_dir = dataset_dir.replace('\\', '/')
    dir_els = dataset_dir.split('/')

    if len(dir_els) < 2:
        print(f'Cannot determine species and dataset from pathname {dataset_dir}')

    dataset = dir_els[-1]
    species = dir_els[-2]

    db_file = os.path.join(dataset_dir, 'db.sqlite3')

    if os.path.isfile(db_file):
        os.remove(db_file)

    engine = create_engine('sqlite:///' + db_file, echo=False, poolclass=NullPool)
    Base.metadata.create_all(engine)
    db_connection = engine.connect()
    engine.session = Session(bind=db_connection)
    session = engine.session

    save_genomic_dataset_details(session, species, dataset)

    if 'Reference_sets' in study_data:
        for file in study_data['Reference_sets']:
            update_genomic_ref(session, os.path.join(dataset_dir, file))
            print(f'Processed reference set {file} for species {species} dataset {dataset}')
    else:
        print('Error - no reference sets were provided for the dataset.')
        return

    needed_study_items = {'Study', 'Date', 'Institute', 'Study_description', 'Researcher', 'Reference', 'Contact'}
    needed_subject_items = {'Name_in_study', 'Age', 'Sex', 'Annotation_file', 'Annotation_format', 'Annotation_method', 'Annotation_reference'}
    needed_assembly_items = {'Assembly_id', 'Assembly_reference', 'Assembly_file', 'Chromosome', 'Start_CoOrd', 'End_CoOrd'}


    for _, study in study_data['Studies'].items():
        if needed_study_items - set(list(study.keys())):
            print(f'Error - study attributes missing: {",".join(list(needed_study_items - set(study.keys())))}')
            return

        study_obj = save_genomic_study(session, study['Study'], study['Date'], study['Institute'], study['Study_description'], study['Researcher'],
                           study['Reference'], study['Contact'])

        for name, subject in study['Subjects'].items():
            if needed_subject_items - set(subject.keys()):
                print(f'Error - subject attributes missing: {",".join(list(needed_subject_items - set(subject.keys())))}')
                return

            report_link = '/'.join(['study_data', 'Genomic', species.replace(' ', '_'), dataset.replace(' ', '_'), subject['Annotation_file']])
            subject_obj = save_genomic_subject(name, subject['Name_in_study'], subject['Age'], subject['Sex'], report_link,
                            subject['Annotation_method'], subject['Annotation_format'], subject['Annotation_reference'], study_obj)

            for _, assembly in subject['Assemblies'].items():
                if needed_assembly_items - set(assembly.keys()):
                    print(f'Error - assembly attributes missing: {",".join(list(needed_assembly_items - set(assembly.keys())))}')
                    return

                if '>' in assembly['Assembly_file']:
                    filename, seqname = assembly['Assembly_file'].split('>')
                    seqs = simple.read_fasta(os.path.join(dataset_dir, filename))

                    if seqname not in seqs:
                        print(f'Sequence {seqname} not found in assembly file {filename}')
                        return

                    assembly_seq = seqs[seqname]
                else:
                    assembly_seq = simple.read_single_fasta(os.path.join(dataset_dir, assembly['Assembly_file']))

                assembly_obj = save_genomic_assembly(assembly['Assembly_id'], assembly['Assembly_reference'],
                                                     assembly['Assembly_file'], assembly_seq, assembly['Chromosome'], assembly['Start_CoOrd'],
                                                     assembly['End_CoOrd'], subject_obj)

                if subject_obj.annotation_format == 'IMGT':
                    # process_imgt_assembly(assembly_obj)
                    pass
                elif subject_obj.annotation_format == 'VDJbase':
                    process_digger_record(session, study_obj, assembly_obj, dataset_dir, subject_obj, subject['Annotation_file'])
                else:
                    print('%s: Invalid type/format %s' % (dataset_file, assembly['Annotation_format']))

    build_gff(session, dataset_dir)
    return

