# Create database and associated files for a single genomic dataset
from datetime import date
import json

from receptor_utils import simple_bio_seq as simple
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, exc
from sqlalchemy.pool import NullPool

import os
import yaml

from db.genomic_ref import update_genomic_ref, read_gene_order
from db.genomic_airr_model import Sample, Study, Patient, SeqProtocol, TissuePro, DataPro
from db.genomic_db import Base, RefSeq
from db.genomic_db_functions import save_genomic_dataset_details, save_genomic_ref_seq, calculate_appearances, calculate_max_cov_sample
from db.igenotyper import process_genomic_record, add_gene_level_features
from db.bed_file import read_bed_files
from db.source_details import db_source_details


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
        commit_id, branch = db_source_details()
        save_genomic_dataset_details(session, species, dataset, commit_id, branch)

        for val in ('Reference_set_version', 'Reference_sets'):
            if val not in study_data:
                raise ImportException(f'Error - {val} is missing from the study metadata.')
            
        read_gene_order(session, dataset_dir)
        for file in study_data['Reference_sets']:
            update_genomic_ref(session, os.path.join(dataset_dir, file))
            print(f'Processed reference set {file} for species {species} dataset {dataset}')

        reference_set_version = study_data['Reference_set_version']

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
        else:
            # create a dummy reference with no coordinates or features
            save_genomic_ref_seq(session, 'dummy', 'dummy', '', 'dummy', 0, 0, '+')
            ref = session.query(RefSeq).filter(RefSeq.name == 'dummy').one_or_none()
            reference_features = {}

        session.commit()
        for study_name, study in study_data['Studies'].items():
            process_study(dataset_dir, reference_features, session, study, study_name, reference_set_version)

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
    ref_seq_name = ref['name']

    # the ref is always stored in forward sense
    # we record the sense in which it is presented. This is used in conjunction with the sense markers in imported_genes.csv
    if ref['sense'] == '-':
        ref_seqs[ref_seq_name] = simple.reverse_complement(ref_seqs[ref_seq_name])

    simple.write_fasta(os.path.join('samples', species + '_' + ref['sequence_file']), ref_seqs)

    save_genomic_ref_seq(session, ref['name'], ref_seqs[ref_seq_name], ref['reference'], ref['chromosome'], ref['start'], ref['end'], ref['sense'])
    reference_features = read_bed_files(ref['locations'], ref['sense'], len(ref_seqs[ref_seq_name]))
    return reference_features


required_fields = {
    'study': {'study_id', 'accession_reference', 'study_title', 'study_description', 'study_type', 'inclusion_exclusion_criteria', 'grants', 'study_contact', 'collected_by', 'lab_name', 'lab_address', 'submitted_by', 'pub_ids', 'keywords_study'},
    'patient': {'subject_id', 'synthetic', 'species', 'sex', 'age_min', 'age_max', 'age_unit', 'age_event', 'ancestry_population', 'ethnicity', 'race', 'strain_name', 'linked_subjects', 'link_type', 'study_group_description', 'disease_diagnosis', 'disease_length', 'disease_stage', 'prior_therapies', 'immunogen', 'intervention', 'medical_history'},
    'sample': {'sample_id', 'sample_type', 'anatomic_site', 'disease_state_sample', 'collection_time_point_relative', 'collection_time_point_relative_unit', 'collection_time_point_reference', 'biomaterial_provider', 'reference_assembly'},
    'tissuepro': {'tissue_processing', 'tissue', 'cell_subset', 'cell_phenotype', 'cell_species', 'single_cell', 'cell_number', 'cells_per_reaction', 'cell_storage', 'cell_quality', 'cell_isolation', 'cell_processing_protocol'},
    'seqprotocol': {'template_class', 'template_quality', 'template_amount', 'template_amount_unit', 'library_generation_method', 'library_generation_protocol', 'library_generation_kit_version', 'complete_sequences', 'physical_linkage', 'pcr_target_locus', 'forward_pcr_primer_target_location', 'reverse_pcr_primer_target_location', 'sequencing_run_id', 'total_reads_passing_qc_filter', 'sequencing_platform', 'sequencing_facility', 'sequencing_run_date', 'sequencing_kit', 'read_length', 'paired_read_length'},
    'datapro': {'sequencing_data_id', 'primary_annotation', 'file_type', 'filename', 'read_direction', 'paired_filename', 'paired_read_direction', 'software_versions', 'paired_reads_assembly', 'quality_thresholds', 'primer_match_cutoffs', 'collapsing_method', 'data_processing_protocols', 'germline_database', 'data_processing_files', 'analysis_provenance_id'},
}


def check_required_fields(fields, row):
    if required_fields[fields] - set(list(row.keys())):
        raise ImportException(f'Error - metadata attributes missing: {",".join(list(required_fields[fields] - set(list(row.keys()))))}')


def process_study(dataset_dir, reference_features, session, study, study_name, reference_set_version):
    study_obj = None
    tissuepros = []
    seqprotocols = []
    datapros = []
    subject_num = 1
    subjects = {}
    subject_samples = {}

    for val in ('annotation_method', 'annotation_reference'):
        if val not in study:
            raise ImportException(f'Error - {val} is missing from the study metadata.')

    annotation_method = study['annotation_method']
    annotation_reference = study['annotation_reference']

    metadata = simple.read_csv(os.path.join(dataset_dir, study['metadata_file']))

    for row in metadata:
        for k, v in row.items():
            if v == 'NA':
                row[k] = ''

        if not study_obj:
            study_obj = create_study(session, study_name, required_fields, row)

        if row['subject_id'] not in subjects:
            if 'vdjbase_name' in row:
                subject_name = '_'.join(row['vdjbase_name'].split('_')[:-1])
            else:
                subject_name = f'{study_name}_I{subject_num}'
                subject_num += 1
                subject_samples[row['subject_id']] = 1

            subject_object = create_subject(session, study_obj, subject_name, row)
            subjects[row['subject_id']] = subject_object
        else:
            subject_object = subjects[row['subject_id']]
            subject_samples[row['subject_id']] += 1

        if 'vdjbase_name' in row:
            sample_name = row['vdjbase_name']
        else:
            sample_name = f'{subject_object.patient_name}_S{subject_samples[row["subject_id"]]}'

        tissuepro_object = find_or_create_tissuepro(session, tissuepros, row)
        seqprotocol_object = find_or_create_seqprotocol(session, seqprotocols, row)
        datapro_object = find_or_create_datapro(session, datapros, row)
        sample_obj = create_sample(session, study_obj, subject_object, sample_name, tissuepro_object, seqprotocol_object, datapro_object, row, annotation_method, annotation_reference)

        process_genomic_record(session, dataset_dir, sample_obj, study['annotation_file'], reference_features, study['bam_dir'])
        
    session.commit()
    calculate_appearances(session)
    calculate_max_cov_sample(session)
    session.commit()


def create_subject(session, study_obj, subject_name, row):
    check_required_fields('patient', row)

    row['synthetic'] = 'T' in row['synthetic'].upper()

    try:
        species_rec = json.loads(row['species'])
        species_id = species_rec['id']
        species_label = species_rec['label']
    except json.JSONDecodeError:
        raise ImportException('Error - species is not a valid ontology object.')

    if species_label and not species_id or ':' not in species_id:
        print(f'Invalid species ID. for {species_label} id {species_id}')

    row['sex'] = row['sex'].upper()
    if row['sex'] in ('F', 'FEMALE'):
        row['sex'] = 'female'
    elif row['sex'] in ('M', 'MALE'):
        row['sex'] = 'male'
    elif row['sex'] == '':
        row['sex'] = 'not collected'
    else:
        raise ImportException(f"study {study_obj.study_name} sample {row['subject_id']}: invalid sex: {row['sex']}")
        
    try:
        age_min = age_max = None
        if row['age_min']:
            age_min = float(row['age_min'])
        if row['age_max']:
            age_max = float(row['age_max'])
    except ValueError:
        raise ImportException('Error - age_min and age_max must be numeric.')
        
    try:
        age_unit_id = age_unit_label = None
        if row['age_unit']:
            age_unit_rec = json.loads(row['age_unit'])
            age_unit_id = age_unit_rec['id']
            age_unit_label = age_unit_rec['label']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - age_unit is not a valid ontology object: {e}')

    if age_unit_label and not age_unit_id or ':' not in age_unit_id:
        print(f'Error - invalid age_unit ID. for {age_unit_label} id {age_unit_id}')

    mother_in_study = None
    father_in_study = None
    if row['linked_subjects']:
        linked_subjects = row['linked_subjects'].split(',')
        link_type = row['link_type'].split(',')

        if len(linked_subjects) != len(link_type):
            raise ImportException('Error - linked_subjects and link_type must have the same number of entries.')
            
        linked_subjects = [x.strip() for x in linked_subjects]
        link_type = [x.strip() for x in link_type]

        for sub, stype in zip(linked_subjects, link_type):
            if stype.lower() == 'mother':
                mother_in_study = sub
            elif stype.lower() == 'father':
                father_in_study = sub

    try:
        disease_diagnosis_id = disease_diagnosis_label = None
        if row['disease_diagnosis']:
            disease_rec = json.loads(row['disease_diagnosis'])
            disease_diagnosis_id = disease_rec['id']
            disease_diagnosis_label = disease_rec['label']
    except json.JSONDecodeError:
        raise ImportException('Error - disease_diagnosis is not a valid ontology object.')

    if disease_diagnosis_label and not disease_diagnosis_id or ':' not in disease_diagnosis_id:
        print(f'Error - invalid disease_diagnosis ID. for species {disease_diagnosis_label} id {disease_diagnosis_id}')

    subject_obj = Patient(
        subject_id=row['subject_id'],
        synthetic=row['synthetic'],
        species_id=species_id,
        species_label=species_label,
        organism_id='',
        organism_label='',
        sex=row['sex'],
        age_min=age_min,
        age_max=age_max,
        age_unit_id=age_unit_id,
        age_unit_label=age_unit_label,
        age_event=row['age_event'],
        ancestry_population=row['ancestry_population'],
        ethnicity=row['ethnicity'],
        race=row['race'],
        strain_name=row['strain_name'],
        linked_subjects=row['linked_subjects'],
        link_type=row['link_type'],
        study_group_description=row['study_group_description'],
        disease_diagnosis_id=disease_diagnosis_id,
        disease_diagnosis_label=disease_diagnosis_label,
        disease_length=row['disease_length'],
        disease_stage=row['disease_stage'],
        prior_therapies=row['prior_therapies'],
        immunogen=row['immunogen'],
        intervention=row['intervention'],
        medical_history=row['medical_history'],
        patient_name=subject_name,
        mother_in_study=mother_in_study,
        father_in_study=father_in_study,
        study_id=study_obj.id,
    )
    session.add(subject_obj)
    session.commit()
    return subject_obj


def create_study(session, study_name, required_fields, row):
    check_required_fields('study', row)

    try:
        type_rec = json.loads(row['study_type'])
        study_type_id = type_rec['id']
        study_type_label = type_rec['label']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - study type is not a valid ontology object: {e}')

    if study_type_label and not study_type_id or ':' not in study_type_id:
        print(f'Error - invalid study_type ID. for {study_type_label} id {study_type_id}')

    try:
        if row['keywords_study']:
            keywords_study = json.loads(row['keywords_study'])
            keywords_study = ', '.join(keywords_study)
        else:
            keywords_study = ''
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - keywords_study is not a valid JSON list: {e}')

    study_obj = Study(
                study_id=row['study_id'],
                study_title=row['study_title'],
                study_type_id=study_type_id,
                study_type_label=study_type_label,
                study_description=row['study_description'],
                inclusion_exclusion_criteria=row['inclusion_exclusion_criteria'],
                grants=row['grants'],
                study_contact=row['study_contact'],
                collected_by=row['collected_by'],
                lab_name=row['lab_name'],
                lab_address=row['lab_address'],
                submitted_by=row['submitted_by'],
                pub_ids=row['pub_ids'],
                keywords_study=keywords_study,
                adc_publish_date='',
                adc_update_date='',
                num_subjects=0,
                num_samples=0,
                accession_reference=row['accession_reference'],
                study_name=study_name,
            )
    session.add(study_obj)
    session.commit()
    return study_obj


def find_or_create_tissuepro(session, tissuepros, row):
    check_required_fields('tissuepro', row)

    tissuepro_dict = {k: v for k, v in row.items() if k in required_fields['tissuepro']}

    try:
        type_rec = json.loads(row['tissue'])
        tissuepro_dict['tissue_id'] = type_rec['id']
        tissuepro_dict['tissue_label'] = type_rec['label']
        del tissuepro_dict['tissue']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - tissue is not a valid ontology object: {e}')

    if tissuepro_dict['tissue_label'] and not tissuepro_dict['tissue_id'] or ':' not in tissuepro_dict['tissue_id']:
        print(f"Error - invalid tissue ID. for tissue {tissuepro_dict['tissue_label']} id {tissuepro_dict['tissue_id']}")

    try:
        tissuepro_dict['cell_species_id'] = tissuepro_dict['cell_species_label'] = None
        if row['cell_species']:
            type_rec = json.loads(row['cell_species'])
            tissuepro_dict['cell_species_id'] = type_rec['id']
            tissuepro_dict['cell_species_label'] = type_rec['label']
        del tissuepro_dict['cell_species']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - cell_species is not a valid ontology object: {e}')

    if tissuepro_dict['cell_species_label'] and not tissuepro_dict['cell_species_id'] or ':' not in tissuepro_dict['cell_species_id']:
        print(f"Error - invalid cell_species ID. for {tissuepro_dict['cell_species_label']} id {tissuepro_dict['cell_species_id']}")

    try:
        tissuepro_dict['cell_subset_id'] = tissuepro_dict['cell_subset_label'] = None
        if row['cell_subset']:
            type_rec = json.loads(row['cell_subset'])
            tissuepro_dict['cell_subset_id'] = type_rec['id']
            tissuepro_dict['cell_subset_label'] = type_rec['label']
        del tissuepro_dict['cell_subset']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - cell_subset is not a valid ontology object: {e}')

    if tissuepro_dict['cell_subset_label'] and not tissuepro_dict['cell_subset_id'] or ':' not in tissuepro_dict['cell_subset_id']:
        print(f"Error - invalid cell_subset ID. for {tissuepro_dict['cell_subset_label']} id {tissuepro_dict['cell_subset_id']}")


    for obj, td in tissuepros:
        if td == tissuepro_dict:
            return obj
        
    tp_label = f"TP{len(tissuepros) + 1}"

    tissuepro_obj = TissuePro(
        tissue_processing=tp_label,
        tissue_label=tissuepro_dict['tissue_label'],
        cell_subset_id=tissuepro_dict['cell_subset_id'],
        cell_subset_label=tissuepro_dict['cell_subset_label'],
        cell_phenotype=tissuepro_dict['cell_phenotype'],
        cell_species_id=tissuepro_dict['cell_species_label'],
        cell_species_label=tissuepro_dict['cell_species_label'],
        single_cell=tissuepro_dict['single_cell'],
        cell_number=tissuepro_dict['cell_number'],
        cells_per_reaction=tissuepro_dict['cells_per_reaction'],
        cell_storage=tissuepro_dict['cell_storage'],
        cell_quality=tissuepro_dict['cell_quality'],
        cell_isolation=tissuepro_dict['cell_isolation'],
        cell_processing_protocol=tissuepro_dict['cell_processing_protocol'],
    )
    session.add(tissuepro_obj)
    session.commit()
    tissuepros.append([tissuepro_obj, tissuepro_dict])
    return tissuepro_obj


def find_or_create_seqprotocol(session, seqprotocols, row):
    check_required_fields('seqprotocol', row)

    seqprotocol_dict = {k: v for k, v in row.items() if k in required_fields['seqprotocol']}

    try:
        seqprotocol_dict['template_amount_unit_id'] = seqprotocol_dict['template_amount_unit_label'] = None
        if row['template_amount_unit']:
            type_rec = json.loads(row['template_amount_unit'])
            seqprotocol_dict['template_amount_unit_id'] = type_rec['id']
            seqprotocol_dict['template_amount_unit_label'] = type_rec['label']
        del seqprotocol_dict['template_amount_unit']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - template_amount_unit is not a valid ontology object: {e}')

    if seqprotocol_dict['template_amount_unit_label'] and not seqprotocol_dict['template_amount_unit_id'] or ':' not in seqprotocol_dict['template_amount_unit_id']:
        print(f"Error - invalid template_amount ID. for {seqprotocol_dict['template_amount_unit_label']} id {seqprotocol_dict['template_amount_unit_id']}")

    for obj, sd in seqprotocols:
        if sd == seqprotocol_dict:
            return obj

    sp_label = f"SP{len(seqprotocols) + 1}"

    seqprotocol_obj = SeqProtocol(
        template_class=seqprotocol_dict['template_class'],
        template_quality=seqprotocol_dict['template_quality'],
        template_amount=seqprotocol_dict['template_amount'],
        template_amount_unit_id=seqprotocol_dict['template_amount_unit_id'],
        template_amount_unit_label=seqprotocol_dict['template_amount_unit_label'],
        library_generation_method=seqprotocol_dict['library_generation_method'],
        library_generation_protocol=seqprotocol_dict['library_generation_protocol'],
        library_generation_kit_version=seqprotocol_dict['library_generation_kit_version'],
        pcr_target_locus=seqprotocol_dict['pcr_target_locus'],
        forward_pcr_primer_target_location=seqprotocol_dict['forward_pcr_primer_target_location'],
        reverse_pcr_primer_target_location=seqprotocol_dict['reverse_pcr_primer_target_location'],
        complete_sequences=seqprotocol_dict['complete_sequences'],
        physical_linkage=seqprotocol_dict['physical_linkage'],
        sequencing_platform=seqprotocol_dict['sequencing_platform'],
        sequencing_facility=seqprotocol_dict['sequencing_facility'],
        sequencing_kit=seqprotocol_dict['sequencing_kit'],
        read_length=seqprotocol_dict['read_length'],
        paired_read_length=seqprotocol_dict['paired_read_length'],
    )

    session.add(seqprotocol_obj)
    session.commit()
    seqprotocols.append([seqprotocol_obj, seqprotocol_dict])
    return seqprotocol_obj


def find_or_create_datapro(session, datapros, row):
    check_required_fields('datapro', row)

    datapro_dict = {k: v for k, v in row.items() if k in required_fields['datapro']}

    for obj, dd in datapros:
        if dd == datapro_dict:
            return obj

    dp_label = f"DP{len(datapros) + 1}"

    datapro_obj = DataPro(
        data_processing_id=dp_label,
        primary_annotation=datapro_dict['primary_annotation'],
        software_versions=datapro_dict['software_versions'],
        paired_reads_assembly=datapro_dict['paired_reads_assembly'],
        quality_thresholds=datapro_dict['quality_thresholds'],
        primer_match_cutoffs=datapro_dict['primer_match_cutoffs'],
        collapsing_method=datapro_dict['collapsing_method'],
        data_processing_protocols=datapro_dict['data_processing_protocols'],
        data_processing_files=datapro_dict['data_processing_files'],
        germline_database=datapro_dict['germline_database'],
        analysis_provenance_id=datapro_dict['analysis_provenance_id'],
    )

    session.add(datapro_obj)
    session.commit()
    datapros.append([datapro_obj, datapro_dict])
    return datapro_obj


def create_sample(session, study, patient, sample_name, tissuepro, seqprotocol, datapro, row, annotation_method, annotation_reference):
    check_required_fields('sample', row)

    try:
        row['collection_time_point_relative_unit_id'] = row['collection_time_point_relative_unit_label'] = None
        if row['collection_time_point_relative_unit']:
            type_rec = json.loads(row['collection_time_point_relative_unit'])
            row['collection_time_point_relative_unit_id'] = type_rec['id']
            row['collection_time_point_relative_unit_label'] = type_rec['label']
        del row['collection_time_point_relative_unit']
    except json.JSONDecodeError as e:
        raise ImportException(f'Error - template_amount_unit is not a valid ontology object: {e}')

    if row['collection_time_point_relative_unit_label'] and not row['collection_time_point_relative_unit_id'] or ':' not in row['collection_time_point_relative_unit_id']:
        print(f"Error - invalid collection_time_point_relative_unit ID. for  {row['collection_time_point_relative_unit_label']} id {row['collection_time_point_relative_unit_id']}")

    ref_id = None
    if row['reference_assembly'] and 'none' not in row['reference_assembly'].lower():
        try:
            ref = session.query(RefSeq).filter(RefSeq.name == row['reference_assembly']).one_or_none()
        except exc.NoResultFound:
            ref = None

        if ref is None:
            raise ImportException(f'Error - reference assembly {row["reference_assembly"]} not found in database.')
        else:
            ref_id = ref.id

    sample_obj = Sample(
        sample_name=sample_name,
        sample_id=row['sample_id'],
        sample_type=row['sample_type'],
        anatomic_site=row['anatomic_site'],
        disease_state_sample=row['disease_state_sample'],
        collection_time_point_relative=row['collection_time_point_relative'],
        collection_time_point_relative_unit_id=row['collection_time_point_relative_unit_id'],
        collection_time_point_relative_unit_label=row['collection_time_point_relative_unit_label'],
        collection_time_point_reference=row['collection_time_point_reference'],
        biomaterial_provider=row['biomaterial_provider'],
        reference_assembly=row['reference_assembly'],
        annotation_path='',
        ref_seq_id=ref_id,
        study_id=study.id,
        patient_id=patient.id,
        tissue_pro_id=tissuepro.id,
        seq_protocol_id=seqprotocol.id,
        data_pro_id=datapro.id,
        annotation_reference=annotation_reference,
        annotation_method=annotation_method,
    )
    session.add(sample_obj)
    session.commit()
    return sample_obj
