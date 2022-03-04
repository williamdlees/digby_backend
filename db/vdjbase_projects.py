# Code to import projects, driven by the contents of projects.yml
import copy
import os
import os.path
from glob import glob
import yaml
import traceback
import json
import sys
import receptor_utils.simple_bio_seq as simple

from db.vdjbase_exceptions import DbCreationError
from db.vdjbase_airr_model import SeqProtocol, Study, TissuePro, Patient, Sample, DataPro
from db.vdjbase_airr_common import read_definition_data




# a table of compound genes used in the pipeline, eg TRBV3-12, and their equiv in VDJbase eg TRBV3-1/2
# the table is not stored in the database ..yet.. because it's only needed during intake
compound_genes = {}


# Read study, subject, sample metadata from definition files into the database
# The key source of metadata is MiAIRR, in json format. Where this is present, the file airr_correspondence.csv
# contains a mapping between the airr file, airr repertoire id, and the VDJbase Pm_In_Sn.
# The file projects.yml (or individual yml files in each sample folder) contains additional
# data added by the processing pieline, e.g. pipeline software versions, read count, etc.
# If no MiAIRR data is found for a sample, we take what we can from the yml files.
#
# The mapping between MiAIRR/YML fields and the SQLAlchemy classes is contained in db/vdjbase_airr_schema_defs.csv.
# This file is generated from the airr schema by parse_airr_schema.py and manually annotated as vdjbase_airr_schema_defs.xls
# before export to csv. Class and filter definitions are created by vdjbase_create_airr_classes.py
# and vdjbase_create_api_filters.py
#
# Metadata records are created in meta_records. This makes it possible to look up the id of a subsidiary record, eg TissueProcessing,
# by doing dict comparisons againts meta_records, un-flattening the airr schema structure. For a record to be
# shared, all fields (including absent fields) must match. All records are committed to the database once processing
# of MiAIRR AND YML is complete, to remove any dependency on when the data is gathered. A key assumption is that the
# database is built from scratch, hence we know in advance what record ids will be allocated


def import_studies(ds_dir, species, dataset, session):
    result = []
    study_data = {}

    study_file = os.path.join(ds_dir, 'projects.yml')
    if os.path.isfile(study_file):
        with open(study_file, 'r') as fi:
            yml_data = yaml.safe_load(fi)
    else:
        # try to consolidate data from sample-level files
        yml_data = consolidate_metadata(ds_dir)

    # we read this directly so that we retain the '.' delimiters for some fields in simple_name
    table_fields = read_definition_data(remove_dots=False)

    # process MiAIRR metadata if present
    yml_full = True

    if os.path.isfile(os.path.join(ds_dir, 'airr_correspondence.csv')):
        airr_corresp = simple.read_csv(os.path.join(ds_dir, 'airr_correspondence.csv'))
        process_airr_metadata(ds_dir, airr_corresp, table_fields, session)
        yml_full = False

    process_yml_metadata(yml_data, table_fields, yml_full)


# Read airr metadata, using info from the correspondence file. Add to meta_records
# Note that more than one repertoire may contribute to a VDJbase sample if records
# have been combined, hence we need to cope with data merging in any field
def process_airr_metadata(ds_dir, airr_corresp, table_fields, session):
    miairr_json = {}

    for rec in airr_corresp:
        meta_records = {
            'Study': {},
            'Patient': {},
            'TissuePro': {},
            'SeqProtocol': {},
            'DataPro': {},
            'GenoDetection': {},
            'Sample': {}
        }

        if rec['airr_file'] not in miairr_json:
            with open(os.path.join(ds_dir, rec['airr_file']), 'r') as fi:
                miairr_json[rec['airr_file']] = json.load(fi)
        rep_ids = rec['airr_repertoire_id'].split(';')

        # specify tables explicitly, so that they get traversed in the desired order
        tables = ['Study', 'TissuePro', 'SeqProtocol', 'DataPro', 'Patient', 'Sample']
        build_metadata(tables, table_fields, rec['vdjbase_name'], miairr_json, rec['airr_file'], rep_ids, meta_records)
        merge_attributes(meta_records, table_fields)
        commit_database(meta_records, rec['vdjbase_name'], session)


# Build metadata 'rows' for a VDJbase sample by walking the schema, looking for desired attributes
def build_metadata(tables, table_fields, vdjbase_name, miairr_json, miairr_file, rep_ids, meta_records):
    desired_attributes = {}
    desired_compounds = {}

    for table in tables:
        desired_attributes[table] = [rec['simple_name'] for rec in table_fields[table].values() if rec['VDJbase table'] == table]
        desired_compounds[table] = [rec['simple_name'].split('.')[-1] for rec in table_fields[table].values() if rec['VDJbase table'] == table and '.' in rec['simple_name']]

    row = {}

    for rep_id in rep_ids:
        repertoire = find_repertoire(miairr_json[miairr_file], rep_id)

        if not repertoire:
            raise DbCreationError(f'Repertoire id {rep_id} does not exist in MiAIRR file {miairr_file}')

        walk_repertoire(repertoire, desired_attributes, desired_compounds, row, ['Repertoire'], meta_records)


# walk the repertoire looking for attributes on our list
def walk_repertoire(obj, desired_attributes, desired_compounds, row, crumb, meta_records):
    if isinstance(obj, dict):
        for k, v in obj.items():
            walk_repertoire(v, desired_attributes, desired_compounds, row, [*crumb, k], meta_records)
    elif isinstance(obj, list):
        for item in obj:
            walk_repertoire(item, desired_attributes, desired_compounds, row, crumb, meta_records)
    else:
        for table, desired_atts in desired_attributes.items():
            for desired_att in desired_atts:
                if crumb[-1] in desired_compounds[table] and len(crumb) > 1:
                    target = f'{crumb[-2]}.{crumb[-1]}'
                else:
                    target = crumb[-1]

                if target == desired_att:
                    collect_attribute(desired_att, obj, meta_records, table)


# Collect an attribute value into a meta_record
def collect_attribute(desired_att, obj, meta_records, table):
    if desired_att not in meta_records[table]:
        meta_records[table][desired_att] = [obj]
    elif obj not in meta_records[table][desired_att]:
        meta_records[table][desired_att].append(obj)


# merge lists of attrbutes as sensibly as we can, converting to type spec in table_defs
def merge_attributes(meta_records, table_fields):
    for table, row_spec in table_fields.items():
        for attr, vals in list(meta_records[table].items()):
            try:
                if vals == [None]:
                    meta_records[table][attr] = None
                elif row_spec[attr]['type'] == 'string':
                    res = []
                    for v in vals:
                        if v and str(v) not in res:
                            res.append(str(v))
                            meta_records[table][attr] = ','.join(res)
                elif row_spec[attr]['type'] == 'integer':
                    meta_records[table][attr] = 0
                    for v in vals:
                        if v:
                            meta_records[table][attr] += int(v)
                elif row_spec[attr]['type'] == 'number':
                    meta_records[table][attr] = 0
                    for v in vals:
                        if v:
                            meta_records[table][attr] += float(v)
                elif row_spec[attr]['type'] == 'boolean':
                    meta_records[table][attr] = 0
                    for v in vals:
                        if v:
                            meta_records[table][attr] = 1
            except:
                print(f"Error combining value {v} in attribute {attr}: cannot coerce to {row_spec[attr]['type']}")

            if '.' in attr:
                meta_records[table][attr.replace('.', '_')] = meta_records[table][attr]
                del meta_records[table][attr]


# Find a repertoire in the json record for a single MiAIRR file
def find_repertoire(miairr_file, rep_id):
    for rep in miairr_file['Repertoire']:
        if rep['repertoire_id'] == rep_id:
            return rep

    return None


# Commit the records for a single repertoire
def commit_database(meta_records, vdjbase_name, session):
    row_ids = {}
    (p_n, i_n, s_n) = vdjbase_name.split('_')
    meta_records['Study']['name'] = p_n
    meta_records['Patient']['name'] = f'{p_n}_{i_n}'
    meta_records['Sample']['name'] = f'{p_n}_{i_n}_{s_n}'
    meta_records['Sample']['sample_group'] = s_n.replace('S', '')

    for table in [Study, TissuePro, SeqProtocol, DataPro]:
        table_name = table.__name__

        db_row = session.query(table).filter_by(**meta_records[table_name]).one_or_none()
        if db_row is None:
            db_row = table(**meta_records[table_name])
            session.add(db_row)
            session.flush()
        row_ids[table.__name__] = db_row.id

    meta_records['Patient']['study_id'] = row_ids['Study']
    patient = Patient(**meta_records['Patient'])
    session.add(patient)
    session.flush()

    meta_records['Sample']['patient_id'] = patient.id
    meta_records['Sample']['seq_protocol_id'] = row_ids['SeqProtocol']
    meta_records['Sample']['study_id'] = row_ids['Study']
    meta_records['Sample']['tissue_pro_id'] = row_ids['TissuePro']
    sample = Sample(**meta_records['Sample'])
    session.add(sample)
    session.commit()



# enumerate dirs and paths under the specified directory
def listdp(dir):
    dirs = [os.path.split(name)[0] for name in glob(os.path.join(dir, '*/'))]
    return zip([os.path.split(name)[1] for name in dirs], dirs)

# Consolidate metadata from yml files in each sample directory
def consolidate_metadata(export_dir):
    metadata = {}
    results = []
    for project_name, project_path in listdp(os.path.join(export_dir, 'samples')):
        for sample_name, sample_path in listdp(project_path):
            metadata_file = os.path.join(sample_path, sample_name + '.yml')
            if not os.path.isfile(metadata_file):
                print('metadata file %s not found: sample will not be included!' % metadata_file)
                results.append('metadata file %s not found: sample will not be included!' % metadata_file)
                continue

            try:
                with open(metadata_file, 'r') as fi:
                    study_data = yaml.safe_load(fi)

                if project_name not in study_data:
                    raise Exception('Project name %s not in sample %s metadata. Sample will not be included.' % (sample_name, project_name))

                if not isinstance(study_data[project_name]['Accession id'], str):
                    raise Exception('Accession id in metadata for sample %s, project %s is not a string. Sample will not be included.' % (sample_name, project_name))

                # fix up some potentially missing fields

                if study_data[project_name]['Samples'][sample_name]['Sequence Protocol Name'] is None:
                    study_data[project_name]['Samples'][sample_name]['Sequence Protocol Name'] = list(study_data[project_name]['Sequence Protocol'].keys())[0]

                if study_data[project_name]['Samples'][sample_name]['Tissue Processing Name'] is None:
                    study_data[project_name]['Samples'][sample_name]['Tissue Processing Name'] = list(study_data[project_name]['Tissue Processing'].keys())[0]

                if study_data[project_name]['Samples'][sample_name]['Reads'] is None:
                    study_data[project_name]['Samples'][sample_name]['Reads'] = 0

                if study_data[project_name]['Samples'][sample_name]['Chain'] is None:
                    study_data[project_name]['Samples'][sample_name]['Chain'] = 'IGH'

                # fix for fields changed with 16 oct release

                for gk, gd in study_data[project_name]['Genotype Detections'].items():
                    if 'Aligner Tool' not in gd and 'Alignment Tool' in gd:
                            study_data[project_name]['Genotype Detections'][gk]['Aligner Tool'] = gd['Alignment Tool']
                            if gd['Alignment reference v'] == gd['Alignment reference d'] and gd['Alignment reference v'] == gd['Alignment reference j']:
                                study_data[project_name]['Genotype Detections'][gk]['Germline Reference'] = gd['Alignment reference v']
                            else:
                                study_data[project_name]['Genotype Detections'][gk]['Germline Reference'] = '(%s (v), %s (d), %s (j))' %\
                                                     (gd['Alignment reference v'], gd['Alignment reference d'], gd['Alignment reference j'])
                    if isinstance(gd['Single Assignment'], str):
                        study_data[project_name]['Genotype Detections'][gk]['Single Assignment'] = ('t' in gd['Single Assignment'] or 'T' in gd['Single Assignment'])

                for sk, sp in study_data[project_name]['Sequence Protocol'].items():
                    if 'Primer.3.location' in sp and 'Primer 3 location' not in sp:
                        study_data[project_name]['Sequence Protocol'][sk]['Primer 3 location'] = sp['Primer.3.location']
                    if 'Primer.5.location' in sp and 'Primer 5 location' not in sp:
                        study_data[project_name]['Sequence Protocol'][sk]['Primer 5 location'] = sp['Primer.5.location']

                for sk, sp in study_data[project_name]['Subjects'].items():
                    if 'Health.Status' in sp and 'Health Status' not in sp:
                        study_data[project_name]['Subjects'][sk]['Health Status'] = sp['Health.Status']
                    if 'Original.name' in sp and 'Original name' not in sp:
                        study_data[project_name]['Subjects'][sk]['Original name'] = sp['Original.name']

                for tk, tp in study_data[project_name]['Tissue Processing'].items():
                    if 'Cell.Type' in tp and 'Cell Type' not in tp:
                        study_data[project_name]['Tissue Processing'][tk]['Cell Type'] = tp['Cell.Type']
                    if 'Sub.Cell.Type' in tp and 'Sub Cell Type' not in tp:
                        study_data[project_name]['Tissue Processing'][tk]['Sub Cell Type'] = tp['Sub.Cell.Type']

                if project_name not in metadata:
                    metadata[project_name] = copy.deepcopy(study_data[project_name])
                else:
                    if sample_name not in study_data[project_name]['Samples']:
                        raise Exception('Wrong sample name not in in sample %s metadata. Sample will not be included.' % sample_name)

                    metadata[project_name]['Samples'][sample_name] = copy.deepcopy(study_data[project_name]['Samples'][sample_name])

                    for section in ('Sequence Protocol', 'Subjects', 'Tissue Processing', 'Genotype Detections'):
                        for item in study_data[project_name][section].keys():
                            if item not in metadata[project_name][section]:
                                metadata[project_name][section][item] = copy.deepcopy(study_data[project_name][section][item])


            except Exception as e:
                print('Exception processing metadata file %s: %s' % (metadata_file, e))

    with open(os.path.join(export_dir, 'consolidated.yml'), 'w') as fo:
        fo.write(yaml.dump(metadata, default_flow_style=False))

    return metadata


