# Create database and associated files for a single genomic dataset
from datetime import date

from receptor_utils import simple_bio_seq as simple
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool

import os
import yaml

from db.genomic_ref import update_genomic_ref, read_gene_order
from db.genomic_db import Base, RefSeq
from db.genomic_db_functions import save_genomic_study, save_genomic_subject, save_genomic_sample, save_genomic_dataset_details, save_genomic_ref_seq, calculate_appearances
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
        engine = create_database(dataset_dir)
        session = engine.session
        save_genomic_dataset_details(session, species, dataset)

        read_gene_order(session, dataset_dir)
        if 'Reference_sets' in study_data:
            for file in study_data['Reference_sets']:
                update_genomic_ref(session, os.path.join(dataset_dir, file))
                print(f'Processed reference set {file} for species {species} dataset {dataset}')
        else:
            raise ImportException('Error - no reference sets were provided for the dataset.')

        reference_set_version = ''
        if 'Reference_set_version' in study_data:
            reference_set_version = study_data['Reference_set_version']
        else:
            raise ImportException('Error - no reference set version was provided for the dataset.')

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
        for study_name, study in study_data['Studies'].items():
            process_study(dataset, dataset_dir, reference_features, session, species, study, study_name, reference_set_version)

    except ImportException as e:
        print(e)
        return

    session.commit()        # final commit in case something is left hanging


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
    return engine


def process_reference_assembly(session, ref, species):
    needed_reference_assembly_items = {'locations', 'sequence_file', 'chromosome', 'start', 'end', 'name', 'reference', 'sense'}

    if needed_reference_assembly_items - set(list(ref.keys())):
        raise ImportException(f'Error - study attributes missing: {",".join(list(needed_reference_assembly_items - set(ref.keys())))}')

    if not os.path.isfile(ref['sequence_file']):
        raise ImportException(f"Error: reference sequence file {ref['sequence_file']} not found.")

    ref_seqs = simple.read_fasta(ref['sequence_file'])
    ref_seq_name = list(ref_seqs.keys())[0]

    # the ref is always stored in forward sense
    # we record the sense in which it is presented. This is used in conjunction with the sense markers in imported_genes.csv
    if ref['sense'] == '-':
        ref_seqs[ref_seq_name] = simple.reverse_complement(ref_seqs[ref_seq_name])

    simple.write_fasta(os.path.join('samples', species + '_' + ref['sequence_file']), ref_seqs)

    save_genomic_ref_seq(session, ref['name'], ref_seqs[ref_seq_name], ref['reference'], ref['chromosome'], ref['start'], ref['end'], ref['sense'])
    reference_features = read_bed_files(ref['locations'], ref['sense'], len(ref_seqs[ref_seq_name]))
    return reference_features


def process_study(dataset, dataset_dir, reference_features, session, species, study, study_name, reference_set_version):
    needed_study_items = {'Study', 'Id', 'Date', 'Institute', 'Study_description', 'Researcher', 'Reference', 'Contact'}
    if needed_study_items - set(list(study.keys())):
        raise ImportException(f'Error - study attributes missing: {",".join(list(needed_study_items - set(study.keys())))}')

    study_date = study['Date']
    if isinstance(study_date, str):
        study_date = date.fromisoformat(study_date)

    study_obj = save_genomic_study(session,
                                   study_name,
                                   study['Study'],
                                   study['Id'],
                                   study_date,
                                   study['Institute'],
                                   study['Study_description'],
                                   study['Researcher'],
                                   study['Reference'],
                                   study['Contact']
                                   )
    session.commit()
    for subject_name, subject in study['Subjects'].items():
        needed_subject_items = {'Name_in_study'}
        if needed_subject_items - set(subject.keys()):
            raise ImportException(f'Error - subject attributes missing: {",".join(list(needed_subject_items - set(subject.keys())))}')

        subject_obj = save_genomic_subject(subject_name, subject['Name_in_study'], study_obj)

        optional_subject_items = [
            ('Age', 'age'),
            ('Sex', 'sex'),
            ('Self-reported Ethnicity', 'self_ethnicity'),
            ('Grouped Ethnicity', 'grouped_ethnicity'),
            ('Population', 'population'),
            ('Pop', 'population_abbr'),
            ('Superpopulation', 'super_population'),
            ('Name_in_study', 'name_in_study'),
            ('Mother_in_study', 'mother_in_study'),
            ('Father_in_study', 'father_in_study'),
        ]

        for yml_attr, sql_attr in optional_subject_items:
            if yml_attr in subject and subject[yml_attr]:
                setattr(subject_obj, sql_attr, subject[yml_attr])

        for name, sample in study['Samples'].items():
            if sample['Subject_name'] == subject_name:
                needed_sample_items = {'Sample_name_in_study', 'Subject_name', 'Reference_assembly', 'Annotation_format'}
                if needed_sample_items - set(sample.keys()):
                    raise ImportException(f'Error - sample attributes missing: {",".join(list(needed_sample_items - set(sample.keys())))}')

                sample_obj = save_genomic_sample(session,
                                                name, 
                                                subject_obj,
                                                sample['Sample_name_in_study'], 
                                                sample['Reference_assembly'], 
                                                sample['Annotation_format'],
                                                ) 

                optional_sample_items = [
                    ('Annotation_file', 'annotation_file'),
                    ('Annotation_method', 'annotation_method'),
                    ('Annotation_reference', 'annotation_reference'),
                    ('Reference_set_version', 'reference_set_version'),
                    ('Locus_coverage', 'locus_coverage'),
                    ('Sequencing_platform', 'sequencing_platform'),
                    ('Assembly_method', 'assembly_method'),
                    ('DNA_source', 'dna_source'),
                ]

                for yml_attr, sql_attr in optional_sample_items:
                    if yml_attr in sample and sample[yml_attr]:
                        setattr(sample_obj, sql_attr, sample[yml_attr])

                # IMGT and Digger formats not currently supported
                if sample_obj.annotation_format == 'IGenotyper':
                    if not reference_features:
                        raise ImportException('Error: Igenotyper records require a reference assembly and one or more bed files')

                    process_igenotyper_record(session, dataset_dir, sample_obj, sample['Annotation_file'], reference_features)
                else:
                    raise ImportException('Error: in yml file: Invalid type/format %s' % (sample_obj.annotation_format))

    session.commit()
    calculate_appearances(session)
    session.commit()
