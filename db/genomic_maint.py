# Create database and associated files for a single genomic dataset
import importlib
import shutil
from datetime import date
from hashlib import sha256

from receptor_utils import simple_bio_seq as simple
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool

import os
import yaml

from db.genomic_ref import update_genomic_ref, read_gene_order
from db.build_gff import build_gff
from db.digger_assembly import process_digger_record
from db.genomic_db import Base, RefSeq
from db.genomic_db_functions import save_genomic_study, save_genomic_subject, save_genomic_assembly, save_genomic_dataset_details, save_genomic_ref_seq, \
    add_feature_to_ref
#from db.imgt_assembly import process_imgt_assembly
from db.igenotyper import process_igenotyper_record, add_gene_level_features
from db.bed_file import read_bed_files

Session = sessionmaker()

class ImportException(Exception):
    pass


def create_dataset(species, dataset):
    try:
        dataset_dir = os.getcwd()

        if not os.path.isdir(os.path.join(dataset_dir, 'samples')):
            os.mkdir(os.path.join(dataset_dir, 'samples'))

        if not os.path.isdir(dataset_dir):
            raise ImportException(f'Directory {dataset_dir} does not exist.')

        study_data = read_yml_file(dataset_dir)
        session = create_database(dataset_dir)
        save_genomic_dataset_details(session, species, dataset)

        read_gene_order(session, dataset_dir)
        if 'Reference_sets' in study_data:
            for file in study_data['Reference_sets']:
                update_genomic_ref(session, os.path.join(dataset_dir, file))
                print(f'Processed reference set {file} for species {species} dataset {dataset}')
        else:
            raise ImportException('Error - no reference sets were provided for the dataset.')

        reference_features = None
        if 'Reference_assemblies' in study_data:
            for ref in study_data['Reference_assemblies'].values():
                reference_features = process_reference_assembly(session, ref, species)

            # once all reference sequences and their bed files are read in, create gene level
            # features for each reference sequence. This allows for the possibility that
            # some bed files contain features for >1 reference sequence

            for ref in study_data['Reference_assemblies'].values():
                ref = session.query(RefSeq).filter(RefSeq.name == ref['name']).one_or_none()
                add_gene_level_features(session, ref, reference_features)

        session.commit()
        for study in study_data['Studies'].values():
            process_study(dataset, dataset_dir, reference_features, session, species, study)

        build_gff(session, dataset_dir)

    except ImportException as e:
        print(e)
        return


def read_yml_file(dataset_dir):
    yml_files = []
    for entry in os.scandir(dataset_dir):
        if entry.is_file() and os.path.splitext(entry.name)[1] in ['.yml', '.yaml']:
            yml_files.append(entry.name)
    if len(yml_files) == 0:
        raise ImportException(f'Error: no yml file found in directory {dataset_dir}.')
    elif len(yml_files) > 1:
        raise ImportException(f'Error: multiple yml files found in directory {dataset_dir}.')
    dataset_file = yml_files[0]
    with open(dataset_file, 'r') as fi:
        study_data = yaml.safe_load(fi)
    return study_data


def create_database(dataset_dir):
    db_file = os.path.join(dataset_dir, 'db.sqlite3')
    if os.path.isfile(db_file):
        os.remove(db_file)
    engine = create_engine('sqlite:///' + db_file, echo=False, poolclass=NullPool)
    Base.metadata.create_all(engine)
    db_connection = engine.connect()
    engine.session = Session(bind=db_connection)
    session = engine.session
    return session


def process_reference_assembly(session, ref, species):
    needed_reference_assembly_items = {'locations', 'sequence_file', 'chromosome', 'start', 'end', 'name', 'reference', 'sense'}

    if needed_reference_assembly_items - set(list(ref.keys())):
        raise ImportException(f'Error - study attributes missing: {",".join(list(needed_reference_assembly_items - set(ref.keys())))}')

    if not os.path.isfile(ref['sequence_file']):
        raise ImportException(f"Error: reference sequence file {ref['sequence_file']} not found.")

    ref_seqs = simple.read_fasta(ref['sequence_file'])
    ref_seq_name = list(ref_seqs.keys())[0]

    if ref['sense'] == '-':
        ref_seqs[ref_seq_name] = simple.reverse_complement(ref_seqs[ref_seq_name])

    simple.write_fasta(ref_seqs, os.path.join('samples', species + '_' + ref['sequence_file']))

    save_genomic_ref_seq(session, ref['name'], ref_seqs[ref_seq_name], ref['reference'], ref['chromosome'], ref['start'], ref['end'])
    reference_features = read_bed_files(ref['locations'], ref['sense'], len(ref_seqs[ref_seq_name]))
    return reference_features


lookup_feature_type = {
    'gene': 'V-GENE',
    'nonamer': 'V-NONAMER',
    'heptamer': 'V-HEPTAMER',
    'spacer': 'V-RSS-SPACER',
}


def process_study(dataset, dataset_dir, reference_features, session, species, study):
    needed_study_items = {'Study', 'Date', 'Institute', 'Study_description', 'Researcher', 'Reference', 'Contact'}
    if needed_study_items - set(list(study.keys())):
        raise ImportException(f'Error - study attributes missing: {",".join(list(needed_study_items - set(study.keys())))}')
    study_obj = save_genomic_study(session, study['Study'], study['Date'], study['Institute'], study['Study_description'],
                                   study['Researcher'],
                                   study['Reference'], study['Contact'])
    session.commit()
    for name, subject in study['Subjects'].items():
        needed_subject_items = {'Name_in_study', 'Age', 'Sex', 'Annotation_file', 'Annotation_format', 'Annotation_method', 'Annotation_reference',
                                'Reference_assembly'}
        if needed_subject_items - set(subject.keys()):
            raise ImportException(f'Error - subject attributes missing: {",".join(list(needed_subject_items - set(subject.keys())))}')

        report_link = '/'.join([species.replace(' ', '_'), dataset.replace(' ', '_'), subject['Annotation_file']])
        subject_obj = save_genomic_subject(session, name, subject['Name_in_study'], subject['Age'], subject['Sex'], report_link,
                                           subject['Annotation_method'], subject['Annotation_format'], subject['Annotation_reference'],
                                           subject['Reference_assembly'], study_obj)

        assembly_objs = []
        for assembly in subject['Assemblies'].values():
            needed_assembly_items = {'Assembly_id', 'Assembly_reference', 'Assembly_file', 'Chromosome', 'Start_CoOrd', 'End_CoOrd'}
            if needed_assembly_items - set(assembly.keys()):
                raise ImportException(f'Error - assembly attributes missing: {",".join(list(needed_assembly_items - set(assembly.keys())))}')

            if '>' in assembly['Assembly_file']:
                filename, seqname = assembly['Assembly_file'].split('>')
                seqs = simple.read_fasta(os.path.join(dataset_dir, filename))

                if seqname not in seqs:
                    raise ImportException(f'Sequence {seqname} not found in assembly file {filename}')

                assembly_seq = seqs[seqname]
            else:
                assembly_seq = simple.read_single_fasta(os.path.join(dataset_dir, assembly['Assembly_file']))

            assembly_obj = save_genomic_assembly(assembly['Assembly_id'], assembly['Assembly_reference'],
                                                 assembly['Assembly_file'], assembly_seq, assembly['Chromosome'], assembly['Start_CoOrd'],
                                                 assembly['End_CoOrd'], subject_obj)

            assembly_objs.append(assembly_obj)

        if subject_obj.annotation_format == 'IMGT':
            # process_imgt_assembly(assembly_obj)
            pass
        elif subject_obj.annotation_format == 'VDJbase':
            for assembly_obj in assembly_objs:
                process_digger_record(session, species, assembly_obj, dataset_dir, subject_obj, subject['Annotation_file'], reference_features)
        elif subject_obj.annotation_format == 'IGenotyper':
            if not reference_features:
                raise ImportException('Error: Igenotyper records require a reference assembly and one or more bed files')

            process_igenotyper_record(session, species, dataset_dir, subject_obj, subject['Annotation_file'], reference_features)
        else:
            raise ImportException('Error: in yml file: Invalid type/format %s' % (assembly['Annotation_format']))


