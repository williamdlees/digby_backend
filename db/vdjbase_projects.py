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
from sqlalchemy import inspect

from db.vdjbase_exceptions import DbCreationError
from db.vdjbase_airr_model import SeqProtocol, Study, TissuePro, Patient, Sample, DataPro, GenoDetection
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
    result = ['Processing metadata']

    study_file = os.path.join(ds_dir, 'projects.yml')
    if os.path.isfile(study_file):
        with open(study_file, 'r') as fi:
            yml_data = yaml.safe_load(fi)
    else:
        # try to consolidate data from sample-level files
        yml_data = consolidate_metadata(ds_dir)

    # we read this directly so that we retain the '.' delimiters for some fields in simple_name
    table_fields = read_definition_data(remove_dots=False)

    miairr_projects = []
    airr_corresp = []

    if os.path.isfile(os.path.join(ds_dir, 'airr_correspondence.csv')):
        airr_corresp = simple.read_csv(os.path.join(ds_dir, 'airr_correspondence.csv'))

        for rec in airr_corresp:
            proj = rec['vdjbase_name'].split('_')[0]
            if proj not in miairr_projects:
                miairr_projects.append(proj)

    for project_name in yml_data.keys():
        miairr_metadata = {}
        if project_name in miairr_projects:
            miairr_metadata = process_airr_metadata(project_name, ds_dir, airr_corresp, table_fields, session)

        process_yml_metadata(project_name, miairr_metadata, yml_data, table_fields, session)

    for project_name in miairr_projects:
        if project_name not in yml_data:
            miairr_metadata = process_airr_metadata(project_name, ds_dir, airr_corresp, table_fields, session)
            process_yml_metadata(project_name, miairr_metadata, yml_data, table_fields, session)

    session.commit()
    return result


def process_yml_metadata(project_name, miairr_metadata, yml_data, table_fields, session):
    sample_names = []
    if project_name in yml_data:
        yml_project_data = yml_data[project_name]
        sample_names = yml_project_data['Samples'].keys()
    else:
        yml_project_data = None
        sample_names = miairr_metadata.keys()

    if not sample_names:
        print(f"Error: No MiAIRR data or YML data found for project {project_name}")
        return

    for sample_name in sample_names:
        meta_records = {
            'Study': {},
            'Patient': {},
            'TissuePro': {},
            'SeqProtocol': {},
            'DataPro': {},
            'GenoDetection': {},
            'Sample': {}
        }

        if sample_name not in miairr_metadata:
            if yml_project_data is None:
                print(f"Error: No MiAIRR data or YML data found for sample {sample_name}")
                continue

            # get everything we can from yml
            build_yml_metadata(sample_name, table_fields, yml_project_data, meta_records, True)
            commit_database(meta_records, sample_name, session)     # fall back to YML exclusively
            continue

        if yml_project_data is not None:
            build_yml_metadata(sample_name, table_fields, yml_project_data, meta_records, False)

        def setattrs(_self, **kwargs):
            for k, v in kwargs.items():
                _self[k] = v

        for table in meta_records.keys():
            if table in miairr_metadata[sample_name]:
                setattrs(meta_records[table], **miairr_metadata[sample_name][table])

        commit_database(meta_records, sample_name, session)

# apply some controls to fields that often seem to go awry
def fixup_fields(meta_records):
    if 'sex' in meta_records['Patient'] and meta_records['Patient']['sex']:
        disc = meta_records['Patient']['sex'].upper()[0]
        if disc == 'F':
            meta_records['Patient']['sex'] = 'female'
        elif disc == 'M':
            meta_records['Patient']['sex'] = 'male'
        else:
            meta_records['Patient']['sex'] = None

    if 'tissue_label' in meta_records['TissuePro'] and meta_records['TissuePro']['tissue_label']:
        meta_records['TissuePro']['tissue_label'] = meta_records['TissuePro']['tissue_label'].lower()

    if 'complete_sequences' in meta_records['SeqProtocol']:
        val = str(meta_records['SeqProtocol']['complete_sequences']).lower()
        if 'full' in val or ('complete' in val and 'incomplete' not in val):
            meta_records['SeqProtocol']['complete_sequences'] = 'complete'
        elif 'partial' in val or 'biomed' in val or 'short' in val:
            meta_records['SeqProtocol']['complete_sequences'] = 'partial'

    if 'template_class' in meta_records['SeqProtocol']:
        val = str(meta_records['SeqProtocol']['template_class']).upper()
        if val not in ['RNA', 'DNA']:
            meta_records['SeqProtocol']['template_class'] = None
        else:
            meta_records['SeqProtocol']['template_class'] = val

    if 'library_generation_method' in meta_records['SeqProtocol']:
        val = str(meta_records['SeqProtocol']['library_generation_method']).lower()
        if val not in enum_fields['library_generation_method']:
            meta_records['SeqProtocol']['library_generation_method'] = 'other'


# check that fields that are supposed to be enums are actually in the enum list

enum_fields = {
    "keywords_study": [
                "contains_ig",
                "contains_tr",
                "contains_paired_chain",
                "contains_schema_rearrangement",
                "contains_schema_clone",
                "contains_schema_cell",
                "contains_schema_receptor",
                ],

    "sex": [
            "male",
            "female",
            "pooled",
            "hermaphrodite",
            "intersex",
            "null",
            ],

    "pcr_target_locus": [
            "IGH",
            "IGI",
            "IGK",
            "IGL",
            "TRA",
            "TRB",
            "TRD",
            "TRG",
            ],

    "template_class": [
            "DNA",
            "RNA",
            ],

    "library_generation_method": [
            "PCR",
            "RT(RHP)+PCR",
            "RT(oligo-dT)+PCR",
            "RT(oligo-dT)+TS+PCR",
            "RT(oligo-dT)+TS(UMI)+PCR",
            "RT(specific)+PCR",
            "RT(specific)+TS+PCR",
            "RT(specific)+TS(UMI)+PCR",
            "RT(specific+UMI)+PCR",
            "RT(specific+UMI)+TS+PCR",
            "RT(specific)+TS",
            "other",
            ],
            
    "read_direction": [
            "forward",
            "reverse",
            "mixed",
            "null",
            ],
            
    "complete_sequences": [
            "partial",
            "complete",
            "complete+untemplated",
            "mixed",
            ],

    "physical_linkage": [
            "none",
            "hetero_head-head",
            "hetero_tail-head",
            "hetero_prelinked",
            ],

    "file_type": [
            "fasta",
            "fastq",
            "null",
            ],

    "paired_read_direction": [
            "forward",
            "reverse",
            "mixed",
            "null",
            ],
}


def check_enums(meta_records):
    for table, fields in meta_records.items():
        for field, value in fields.items():
            if field == 'keywords_study' and value:
                for kw in value.split(','):
                    if kw not in enum_fields['keywords_study']:
                        print(f"Warning: keyword {kw} in {table} not in enum list")
            elif field in enum_fields and value and value not in enum_fields[field]:
                print(f"Warning: {field} in {table} has value {value} not in enum list")


# Commit data for one sample, creating other associated rows as needed
def build_yml_metadata(sample_name, table_fields, project_data, meta_records, yml_full):
    desired_attributes = {}

    for table in meta_records.keys():
            desired_attributes[table] = {rec['YML attribute']: rec['simple_name'].replace('.', '_')
                                         for rec in table_fields[table].values()
                                         if rec['VDJbase table'] == table and rec['YML attribute'] and (yml_full or not rec['structured_name'])}

    # apply fixups to GenoDetection fields

    for attr, key in desired_attributes['GenoDetection'].items():
        if key.replace('_', '.') in geno_specials:
            desired_attributes['GenoDetection'][attr] = geno_specials[key.replace('_', '.')]

    # dict attributes are not the same as names - so find indexes into the various YML objects

    yml_objects = {
        'Genotype Detections': {},
        'Sequence Protocol': {},
        'Tissue Processing': {},
    }

    for obj_type in yml_objects.keys():
        for obj_val in project_data[obj_type].values():
            yml_objects[obj_type][obj_val['Name']] = obj_val


    for table in meta_records.keys():
        for desired_attribute, table_attribute in desired_attributes[table].items():
            value = None
            try:
                if '.' not in desired_attribute:
                    value = project_data[desired_attribute]

                else:
                    section, desired_attribute = desired_attribute.split('.')

                    if section == 'Genotype Detections':
                        section_name = project_data['Samples'][sample_name]['Genotype Detection Name']
                        named_section = yml_objects[section][section_name]
                        if desired_attribute in named_section:
                            value = named_section[desired_attribute]
                        elif desired_attribute == 'Alignment Tool' and 'Aligner Tool' in named_section:
                            value = named_section['Aligner Tool']       # fudge for old-style naming used in TCRB
                    elif section == 'Samples':
                        value = project_data['Samples'][sample_name][desired_attribute]
                    elif section == 'Sequence Protocol':
                        section_name = project_data['Samples'][sample_name]['Sequence Protocol Name']
                        named_section = yml_objects[section][section_name]
                        value = named_section[desired_attribute]
                    elif section == 'Subjects':
                        value = project_data['Subjects'][project_data['Samples'][sample_name]['Subject Name']][desired_attribute]
                    elif section == 'Tissue Processing':
                        section_name = project_data['Samples'][sample_name]['Tissue Processing Name']
                        named_section = yml_objects[section][section_name]
                        value = named_section[desired_attribute]
                    else:
                        print(f'Unrecognised section name {section} for attribute {desired_attribute}')
            except Exception as e:
                print(f'Parsing error in section name {section} for attribute {desired_attribute}: {e}')

            if value and table_attribute:
                meta_records[table][table_attribute] = value


# Read airr metadata, using info from the correspondence file. Add to meta_records
# Note that more than one repertoire may contribute to a VDJbase sample if records
# have been combined, hence we need to cope with data merging in any field
def process_airr_metadata(project_name, ds_dir, airr_corresp, table_fields, session):
    miairr_json = {}
    project_meta_records = {}

    for rec in airr_corresp:
        if project_name + '_' in rec['vdjbase_name'] and rec['airr_repertoire_id']:
            meta_records = {
                'Study': {},
                'Patient': {},
                'TissuePro': {},
                'SeqProtocol': {},
                'DataPro': {},
                'GenoDetection': {},
                'Sample': {}
            }

            print(f'Processing airr sample {rec["vdjbase_name"]}')

            if rec['airr_file'] not in miairr_json:
                with open(os.path.join(ds_dir, rec['airr_file']), 'r') as fi:
                    miairr_json[rec['airr_file']] = json.load(fi)
            rep_ids = rec['airr_repertoire_id'].split(';')

            # specify tables explicitly, so that they get traversed in the desired order
            tables = ['Study', 'TissuePro', 'SeqProtocol', 'DataPro', 'Patient', 'Sample', 'GenoDetection']
            build_metadata(tables, table_fields, rec['vdjbase_name'], miairr_json, rec['airr_file'], rep_ids, meta_records)
            merge_attributes(meta_records, table_fields)
            project_meta_records[rec['vdjbase_name']] = meta_records

    if not project_meta_records:
        print(f"No MiAIRR data found for repertoire_id {rec['airr_repertoire_id']} in project {project_name}")

    return project_meta_records


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

        # special for preprocessing tool versions
        if len(crumb) > 2 and crumb[-3] == 'preprocessing' and crumb[-2] == 'software_versions':
            collect_attribute('prepro_tool', f'{crumb[-1]}: {obj}', meta_records, 'GenoDetection')


# Collect an attribute value into a meta_record
def collect_attribute(desired_att, obj, meta_records, table):
    if desired_att not in meta_records[table]:
        meta_records[table][desired_att] = [obj]
    elif obj not in meta_records[table][desired_att]:
        meta_records[table][desired_att].append(obj)


# some GenoDetection specials
geno_specials = {
    'Single Assignment': 'single_assignment',
    'aligner.version': 'aligner_ver',
    'Genotyper.Tool': 'geno_tool',
    'Genotyper.Version': 'geno_ver',
    'Haplotyper.Tool': 'haplotype_tool',
    'Haplotyper.Version': 'haplotype_ver',
}


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
                print(f"Error combining value in attribute {attr}: cannot coerce to {row_spec[attr]['type']}")

            if attr in geno_specials:
                meta_records[table][geno_specials[attr]] = meta_records[table][attr]
                del meta_records[table][attr]
            elif '.' in attr:
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
    fixup_fields(meta_records)
    check_enums(meta_records)

    row_ids = {}
    (p_n, i_n, s_n) = vdjbase_name.split('_')
    meta_records['Study']['study_name'] = p_n
    meta_records['Patient']['patient_name'] = f'{p_n}_{i_n}'
    meta_records['Sample']['sample_name'] = f'{p_n}_{i_n}_{s_n}'
    meta_records['Sample']['sample_group'] = s_n.replace('S', '')

    defaults = {}
    for table in [Study, TissuePro, SeqProtocol, GenoDetection, DataPro]:
        defaults[table.__name__] = {}
        for column_name, column in inspect(table).column_attrs.items():
            if column.expression.nullable:
                defaults[table.__name__][column_name] = None

    for table in [Study, TissuePro, SeqProtocol, GenoDetection, DataPro]:
        table_name = table.__name__

        if not meta_records[table_name]:
            db_row = session.query(table).filter_by(**defaults[table_name]).one_or_none()

            if not db_row:
                print(f'Creating blank record for table {table} in sample {vdjbase_name}')
                db_row = table()
                session.add(db_row)
                session.flush()
            row_ids[table.__name__] = db_row.id
        else:
            # add defaults for any unspecified attributes

            for k, v in defaults[table_name].items():
                if k not in meta_records[table_name]:
                    meta_records[table_name][k] = v

            db_row = session.query(table).filter_by(**meta_records[table_name]).all()
            if db_row and len(db_row) > 1:
                breakpoint()

            # something of a corner case, but we have two projects that are missing some samples in the iReceptor+ metadata.
            # Make sure we only create one project record per study, regardless of metadata differences

            if table == Study:
                db_row = session.query(table).filter(Study.study_name == meta_records[table_name]['study_name']).one_or_none()
            else:
                db_row = session.query(table).filter_by(**meta_records[table_name]).one_or_none()

            if db_row is None:
                db_row = table(**meta_records[table_name])
                session.add(db_row)
                session.commit()
            row_ids[table.__name__] = db_row.id

    patient = session.query(Patient).filter(Patient.patient_name == meta_records['Patient']['patient_name']).one_or_none()

    if not patient:
        meta_records['Patient']['study_id'] = row_ids['Study']
        patient = Patient(**meta_records['Patient'])
        session.add(patient)
        session.flush()

    meta_records['Sample']['patient_id'] = patient.id
    meta_records['Sample']['seq_protocol_id'] = row_ids['SeqProtocol']
    meta_records['Sample']['study_id'] = row_ids['Study']
    meta_records['Sample']['tissue_pro_id'] = row_ids['TissuePro']
    meta_records['Sample']['data_pro_id'] = row_ids['DataPro']
    meta_records['Sample']['geno_detection_id'] = row_ids['GenoDetection']
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
                # print('metadata file %s not found: sample will not be included!' % metadata_file)
                results.append('metadata file %s not found: using json data only' % metadata_file)
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


