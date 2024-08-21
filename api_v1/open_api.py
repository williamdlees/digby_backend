import json
import decimal
import os
import datetime as dt
from app import app
from datetime import datetime
from flask import Blueprint, request, jsonify, Response, send_from_directory
from api.vdjbase.vdjbase import get_vdjbase_species, find_datasets, vdjbase_dbs, sample_info_filters, find_vdjbase_samples, rep_sample_bool_values, VDJBASE_SAMPLE_PATH
from api.genomic.genomic import get_genomic_species, get_genomic_datasets, find_genomic_samples, ceil, genomic_sample_filters, get_genomic_db
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype
from db.vdjbase_airr_model import Sample as Airr_Sample
from schema.models import *
from db.genomic_airr_model import Patient, SeqProtocol, TissuePro, DataPro
from db.vdjbase_airr_model import GenoDetection as Airr_GenoDetection, SeqProtocol as Airr_SeqProtocol, Study as Airr_Study, TissuePro as Airr_TissuePro, Patient as Airr_Patient, DataPro as Airr_DataPro
from db.genomic_airr_model import Study as STD
from db.genomic_airr_model import Sample as SAMPLE
from sqlalchemy import func
from os.path import isfile
from pydantic.fields import FieldInfo
from pydantic import BaseModel, ValidationError
from typing import Any, Optional, Union, get_args, get_origin, List, get_type_hints


api_bp = Blueprint('api_v1', __name__)

def custom_jsonify(obj):
    """Custom JSON encoder for special object types."""

    def encode_obj(o):
        if isinstance(o, Enum):
            return o.value
        if isinstance(o, BaseModel):
            return o.dict()
        if isinstance(o, date):
            return o.isoformat()
        if isinstance(o, dict):
            return {k: encode_obj(v) for k, v in o.items()}
        if isinstance(o, list):
            return [encode_obj(i) for i in o]
        return o
    
    return Response(
        json.dumps(encode_obj(obj)),
        mimetype='application/json'
    )


"""Get species list based on type."""

from api.genomic.genomic import SpeciesApi
from flask_restx import Resource

@api_bp.route('/<type>/species', methods=['GET'])
def get_species(type):
    species_api = SpeciesApi(Resource)
    species = species_api.get()
    
    """Get species list based on type."""
    sp = None
    if type == "genomic":
        sp = get_genomic_species()
    else:
        sp = get_vdjbase_species()

    species_list = list(set(sp))
    species_response_obj = []
    for item in species_list:
        ontology_obj = Ontology(label=item)
        species_response_obj.append(ontology_obj)

    species_response_obj = SpeciesResponse(species=species_response_obj)
    try:
        return jsonify(species_response_obj.dict()), 200

    except Exception as e:
        error_response = ErrorResponse(message=str(e))
        return jsonify(error_response.dict()), 500


@api_bp.route('/<type>/datasets/<species>', methods=['GET'])
def get_species_datasets(type, species):
    """Get datasets for a species based on type."""
    data_sets = None
    if type == "genomic":
        data_sets = get_genomic_datasets(species)

    elif type == "airrseq":
        data_sets = find_datasets(species)

    else:
        error_response = ErrorResponse(message="type not exists")
        return jsonify(error_response.dict()), 500

    data_set_list = []
    for data_set in data_sets:
        data_set_obj = Dataset(dataset=data_set["dataset"], locus=data_set["locus"] if type == "genomic" else None, type=type)
        data_set_list.append(data_set_obj)

    data_set_response = DatasetsResponse(datasets=data_set_list)
    try:
        return jsonify(data_set_response.dict()), 200

    except Exception as e:
        error_response = ErrorResponse(message=str(e))
        return jsonify(error_response.dict()), 500


@api_bp.route('/<type>/subjects/<species>/<dataset>', methods=['GET'])
def get_subject_datasets(type, species, dataset):
    """Get subject datasets for a species and dataset based on type."""
    if type == "genomic":
        subjects_list, is_currect = get_genomic_list_subjects(species, dataset)
        if not is_currect:
            error_response = ErrorResponse(message=str(subjects_list))
            return jsonify(error_response.dict()), 500
        
        else:
            dataset_list = []
            for sample in subjects_list.get('samples'):
                subject_identifier = sample['sample_name'].split('_')[0:1]
                subject_identifier = '_'.join(sample['sample_name'].rsplit('_', 1)[:-1])

                subject_dataset_obj = SubjectDataset(id=sample['sample_id'], 
                                                     study_name=sample['study_name'], 
                                                     subject_identifier=subject_identifier, 
                                                     sample_identifier=sample['sample_name'],
                                                     dataset=sample['dataset'])
                dataset_list.append(subject_dataset_obj)
            subject_dataset_response_obj = SubjectDatasetResponse(subject_datasets=dataset_list)
          
            try:
                return jsonify(subject_dataset_response_obj.dict()), 200

            except Exception as e:
                error_response = ErrorResponse(message=str(e))
                return jsonify(error_response.dict()), 500

    elif type == "airrseq":
        subjects_list, is_currect = get_airrseq_list_subjects(species, dataset)
        if not is_currect:
            error_response = ErrorResponse(message=str(subjects_list))
            return jsonify(error_response.dict()), 500
        
        else:
            dataset_list = []
            for sample in subjects_list.get('samples'):
                subject_identifier = sample['sample_name'].split('_')[0:1]
                subject_identifier = '_'.join(sample['sample_name'].rsplit('_', 1)[:-1])

                subject_dataset_obj = SubjectDataset(id=sample['sample_name'], 
                                                     study_name=sample['sample_name'].split('_')[0], 
                                                     subject_identifier=sample['patient_name'], 
                                                     sample_identifier=sample['sample_name'],
                                                     dataset=sample['dataset'])
                dataset_list.append(subject_dataset_obj)
            subject_dataset_response_obj = SubjectDatasetResponse(subject_datasets=dataset_list)
          
            try:
                return jsonify(subject_dataset_response_obj.dict()), 200

            except Exception as e:
                error_response = ErrorResponse(message=str(e))
                return jsonify(error_response.dict()), 500


def get_genomic_list_subjects(species, genomic_datasets):
    """Get a list of genomic subjects for a species and datasets."""
    args = {
        "page_number": 0,
        "page_size": 100,
        "filter": None,
        "sort_by": None,
        "cols": '["sample_name","study_id"]'
    }

    required_cols = json.loads(args['cols']) if 'cols' in args and args['cols'] else list(genomic_sample_filters.keys())

    for col in required_cols:
        if col not in genomic_sample_filters.keys():
            return 'Bad filter string %s' % args['filter'], False

    if 'study_name' not in required_cols:
        required_cols = ['study_name'] + required_cols
    if 'dataset' not in required_cols:
        required_cols.append('dataset')
    if 'sample_id' not in required_cols:
        required_cols.append('sample_id')
    if 'annotation_path' in required_cols:
        if 'annotation_reference' not in required_cols:
            required_cols.append('annotation_reference')
        if 'annotation_method' not in required_cols:
            required_cols.append('annotation_method')
        if 'contig_bam_path' not in required_cols:
            required_cols.append('contig_bam_path')

    attribute_query = [genomic_sample_filters['sample_id']['field']]        # the query requires the first field to be from Sample

    for col in required_cols:
        if col != 'sample_id' and genomic_sample_filters[col]['field'] is not None:
            attribute_query.append(genomic_sample_filters[col]['field'])

    filter = json.loads(args['filter']) if args['filter'] else []
    datasets = genomic_datasets.split(',')
    ret = find_genomic_samples(attribute_query, species, datasets, filter)

    uniques = {}
    for f in required_cols:
        uniques[f] = []
    uniques['dataset'] = datasets

    # special column for names by dataset

    uniques['names_by_dataset'] = {}
    filter_applied = len(filter) > 0
    if filter_applied:
        for dataset in uniques['dataset']:
            uniques['names_by_dataset'][dataset] = []

    def num_sort_key(x):
        if x is None or x == '':
            return -1
        else:
            try:
                return float(x)
            except:
                return 0

    def name_sort_key(name):
        name = name.split('_')
        for i in range(len(name)):
            name[i] = name[i][1:].zfill(4)
        return name

    for s in ret:
        for f in required_cols:
            if genomic_sample_filters[f]['field'] is not None and 'no_uniques' not in genomic_sample_filters[f]:
                el = s[f]
                if isinstance(el, datetime):
                    el = el.date().isoformat()
                elif isinstance(el, str) and len(el) == 0:
                    el = '(blank)'
                if el not in uniques[f]:
                    uniques[f].append(el)
        if filter_applied:
            uniques['names_by_dataset'][s['dataset']].append(s['sample_id'])

    for f in required_cols:
        try:
            if 'sort' in genomic_sample_filters[f] and genomic_sample_filters[f]['sort'] == 'numeric':
                uniques[f].sort(key=num_sort_key)
            elif 'sort' in genomic_sample_filters[f] and genomic_sample_filters[f]['sort'] == 'underscore':
                uniques[f].sort(key=name_sort_key)
            else:
                uniques[f].sort(key=lambda x: (x is None or x == '', x))
        except:
            pass

    sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None) else []
    if len(sort_specs) == 0:
        sort_specs = [{'field': 'sample_identifier', 'order': 'asc'}]

    for spec in sort_specs:
        f = spec['field']
        if f in genomic_sample_filters.keys():
            if 'sort' in genomic_sample_filters[f] and genomic_sample_filters[f]['sort'] == 'underscore':
                ret = sorted(ret, key=lambda x: name_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
            elif 'sort' in genomic_sample_filters[f] and genomic_sample_filters[f]['sort'] == 'numeric':
                ret = sorted(ret, key=lambda x: num_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
            else:
                ret = sorted(ret, key=lambda x: ((x[f] is None or x[f] == ''),  x[f]), reverse=(spec['order'] == 'desc'))

    total_size = len(ret)

    if args['page_size']:
        first = (args['page_number']) * args['page_size']
        ret = ret[first : first + args['page_size']]

    return {
        'samples': ret,
        'uniques': uniques,
        'total_items': total_size,
        'page_size': args['page_size'],
        'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
    }, True


def get_airrseq_list_subjects(species, dataset):
    """Get a list of AIRR-seq subjects for a species and dataset."""
    if species not in vdjbase_dbs or set(dataset.split(',')).difference(set(vdjbase_dbs[species])):
        return list()

    args = {
        "page_number": 0,
        "page_size": 100,
        "filter": None,
        "sort_by": None,
        "cols": '["sample_name","patient_name"]'
    }

    required_cols = json.loads(args['cols']) if args['cols'] else list(sample_info_filters.keys())
    if 'sample_name' not in required_cols:
        required_cols.append('sample_name')
    if 'dataset' not in required_cols:
        required_cols.append('dataset')

    for col in required_cols:
        if col not in sample_info_filters.keys():
            return 'Bad filter string %s' % args['filter'], False

    if 'genotype' in required_cols:                # needed to compose paths to files
        for field in ('patient_name', 'study_name'):
            if field not in required_cols:
                required_cols.append(field)
        required_cols.append('genotypes')
        required_cols.append('genotype_stats')
        required_cols.append('genotype_report')

    attribute_query = [sample_info_filters['sample_id']['field']]        # the query requires the first field to be from Sample

    for col in required_cols:
        if col != 'sample_id' and sample_info_filters[col]['field'] is not None:
            attribute_query.append(sample_info_filters[col]['field'])

    attribute_query.append(Airr_Sample.id)

    filter = json.loads(args['filter']) if args['filter'] else []
    ret = find_vdjbase_samples(attribute_query, species, dataset.split(','), filter)

    total_size = len(ret)

    uniques = {}

    for f in required_cols:
        uniques[f] = []

    uniques['dataset'] = dataset.split(',')

    # special column for names by dataset

    uniques['names_by_dataset'] = {}
    filter_applied = len(filter) > 0
    if filter_applied:
        for dataset in uniques['dataset']:
            uniques['names_by_dataset'][dataset] = []

    for s in ret:
        for f in required_cols:
            if sample_info_filters[f]['field'] is not None and 'no_uniques' not in sample_info_filters[f]:
                el = s[f]
                if isinstance(el, (dt.datetime, dt.date)):
                    el = el.date().isoformat()
                elif isinstance(el, bool):
                    if f in rep_sample_bool_values:
                        el = rep_sample_bool_values[f][0 if el else 1]
                if (not isinstance(el, str) or len(el) > 0) and el not in uniques[f]:
                    uniques[f].append(el)
                if (sample_info_filters[f]['field'].type.python_type is str) and len(el) == 0 and '(blank)' not in uniques[f]:
                    uniques[f].append('(blank)')

        if filter_applied:
            uniques['names_by_dataset'][s['dataset']].append(s['sample_name'])

    def name_sort_key(name):
        name = name.split('_')
        for i in range(len(name)):
            name[i] = name[i][1:].zfill(4)
        return name

    def num_sort_key(x):
        if x is None or x == '':
            return -1
        else:
            try:
                return float(x)
            except:
                return 0

    for f in required_cols:
        try:
            if 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'underscore':
                uniques[f].sort(key=name_sort_key)
            elif 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'numeric':
                uniques[f].sort(key=num_sort_key)
            elif f != 'names_by_dataset':
                uniques[f].sort(key=lambda x: (x is None or x == '', x))
        except:
            pass

    if 'haplotypes' in required_cols:
        uniques['haplotypes'] = []
        for dset in dataset.split(','):
            session = vdjbase_dbs[species][dset].session
            haplotypes = session.query(HaplotypesFile.by_gene_s).distinct().order_by(HaplotypesFile.by_gene_s).all()
            x = [(h[0]) for h in haplotypes]
            uniques['haplotypes'].extend(x)
        uniques['haplotypes'] = list(set(uniques['haplotypes']))

    sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None) else []
    if len(sort_specs) == 0:
        sort_specs = [{'field': 'name', 'order': 'asc'}]

    for spec in sort_specs:
        f = spec['field']
        if f in sample_info_filters.keys():
            if 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'underscore':
                ret = sorted(ret, key=lambda x: name_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
            elif 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'numeric':
                ret = sorted(ret, key=lambda x: num_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
            else:
                ret = sorted(ret, key=lambda x: ((x[f] is None or x[f] == ''), x[f]), reverse=(spec['order'] == 'desc'))

    if args['page_size']:
        first = (args['page_number']) * args['page_size']
        ret = ret[first:first + args['page_size']]

    for rec in ret:
        for k, v in rec.items():
            if isinstance(v, (dt.datetime, dt.date)):
                rec[k] = v.date().isoformat()
            elif isinstance(v, decimal.Decimal):
                rec[k] = '%0.2f' % v

    if 'genotypes' in required_cols:
        for r in ret:
            r['genotypes'] = {}
            r['genotypes']['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['sample_name'], 'sort_order': 'Locus'})

            r['genotypes']['path'] = app.config['BACKEND_LINK']
            sp = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset']]) + '/'
            r['genotypes']['tigger'] = sp + r['genotype'].replace('samples', '') if r['genotype_stats'] else ''
            r['genotypes']['ogrdbstats'] = sp + r['genotype_stats'].replace('samples', '') if r['genotype_stats'] else ''
            r['genotypes']['ogrdbplot'] = sp + r['genotype_report'].replace('samples', '') if r['genotype_report'] else ''
            del r['genotype_stats']
            del r['genotype_report']

            session = vdjbase_dbs[species][r['dataset']].session
            igsnper_path = session.query(Airr_Sample.igsnper_plot_path, Airr_Sample.genotype_report).filter(Airr_Sample.sample_name == r['sample_name']).one_or_none()

            if igsnper_path is not None and igsnper_path[0] is not None:
                r['genotypes']['igsnper'] = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset'], igsnper_path[0]])
            else:
                r['genotypes']['igsnper'] = ''

    if 'haplotypes' in required_cols:
        for r in ret:
            session = vdjbase_dbs[species][r['dataset']].session
            haplotypes = session.query(Airr_Sample.sample_name, func.group_concat(HaplotypesFile.by_gene_s), func.group_concat(HaplotypesFile.file))
            h = haplotypes.filter(Airr_Sample.sample_name == r['sample_name'])\
                .join(SamplesHaplotype, SamplesHaplotype.samples_id == Airr_Sample.id)\
                .join(HaplotypesFile, HaplotypesFile.id == SamplesHaplotype.haplotypes_file_id)\
                .one_or_none()
            if h is not None and h[1] is not None:
                r['haplotypes'] = {}
                r['haplotypes']['path'] = app.config['BACKEND_LINK']
                for (hap, filename) in zip(h[1].split(','), h[2].split(',')):
                    filename = filename.replace('samples/', '')
                    fp = os.path.join(VDJBASE_SAMPLE_PATH, species, r['dataset'], filename)
                    sp = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset'], filename])
                    if isfile(fp):
                        r['haplotypes'][hap] = {}
                        r['haplotypes'][hap]['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['sample_name'], 'hap_gene': hap, 'sort_order' : 'Locus'})
                        r['haplotypes'][hap]['rabhit'] = sp
            else:
                r['haplotypes'] = ''

    return {
        'samples': ret,
        'uniques': uniques,
        'total_items': total_size,
        'page_size': args['page_size'],
        'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
    }, True

@api_bp.route('/<type>/sample_genotype/<species>/<dataset>/<subject>/<sample>', methods=['GET'])
def get_sample_genotype(type, species, dataset, subject, sample):
    genotype_object = Genotype(receptor_genotype_id='',
                                   locus=Locus('IGH'),
                                   documented_alleles=None,
                                   undocumented_alleles=None,
                                   deleted_genes=None,
                                   inference_process=None)
    
    if type == "genomic":
        return custom_jsonify(genotype_object.dict()), 200

    elif type == "airrseq":
        return custom_jsonify(genotype_object.dict()), 200

    else:
        error_response = ErrorResponse(message=str("type not  exists"))
        return jsonify(error_response.dict()), 500
    

@api_bp.route('/<type>/sample_metadata/<species>/<dataset>/<subject>/<sample>', methods=['GET'])
def get_sample_metadata(type, species, dataset, subject, sample):
    """Get metadata for a specific sample."""
    if type == "genomic":
        subject_info, is_currect = get_genomic_subjet_info(species, dataset, sample)
        if is_currect:
            rep_obj = SampleMetadataResponse(Repertoire=create_repertoire_obj(subject_info))
            return custom_jsonify(rep_obj.dict()), 200

        else:
            error_response = ErrorResponse(message=str(subject_info))
            return jsonify(error_response), 500
        
    elif type == "airrseq":
        subject_info, is_currect = get_airr_subjet_info(species, dataset, sample)
        if is_currect:
            rep_obj = SampleMetadataResponse(Repertoire=create_repertoire_obj(subject_info))
            return custom_jsonify(rep_obj.dict()), 200

        else:
            error_response = ErrorResponse(message=str(subject_info))
            return jsonify(error_response), 500
    else:
        error_response = ErrorResponse(message=str("type not  exists"))
        return jsonify(error_response.dict()), 500


def create_repertoire_obj(subject_info):
    """Create a Repertoire object from subject information."""
    subject_info = fill_missing_required_fields(Repertoire,  subject_info)

    rep_object = Repertoire(repertoire_id=subject_info.get("repertoire_id", None),
                            repertoire_name=subject_info.get("repertoire_name", None),
                            repertoire_description=subject_info.get("repertoire_description", None),
                            study=create_study_object(subject_info),
                            subject=create_subject_objects(subject_info),
                            sample=create_sample_processing_list(subject_info),
                            data_processing=create_data_processing_list(subject_info))
    
    return rep_object


def create_data_processing_list(subject_info):
    """Create a list of DataProcessing objects from subject information."""
    subject_info = fill_missing_required_fields(DataProcessing,  subject_info)

    data_processing_list = []
    data_processing_obj = DataProcessing(data_processing_id=subject_info.get("data_processing_id"),
                                         primary_annotation=subject_info.get("primary_annotation"),
                                         software_versions=subject_info.get("software_versions"),
                                         paired_reads_assembly=subject_info.get("paired_reads_assembly"),
                                         quality_thresholds=subject_info.get("quality_thresholds"),
                                         primer_match_cutoffs=subject_info.get("primer_match_cutoffs"),
                                         collapsing_method=subject_info.get("collapsing_method"),
                                         data_processing_protocols=subject_info.get("data_processing_protocols"),
                                         data_processing_files=[subject_info.get("data_processing_files")] if subject_info.get("data_processing_files") is not None else None,
                                         germline_database=subject_info.get("germline_database"),
                                         germline_set_ref=subject_info.get("germline_set_ref"),
                                         analysis_provenance_id=subject_info.get("analysis_provenance_id"))
    
    data_processing_list.append(data_processing_obj)
    return data_processing_list
    

def create_sample_processing_list(subject_info):
    """Create a list of SampleProcessing objects from subject information."""
    subject_info = fill_missing_required_fields(Sample,  subject_info)
    subject_info = fill_missing_required_fields(CellProcessing,  subject_info)
    subject_info = fill_missing_required_fields(NucleicAcidProcessing,  subject_info)
    subject_info = fill_missing_required_fields(SequencingRun,  subject_info)

    try:
        library_generation_method = LibraryGenerationMethod(subject_info.get("library_generation_method"))
    except:
        library_generation_method = LibraryGenerationMethod("other")

    sample_processing_list = []
    sample_processing_obj = SampleProcessing(sample_processing_id=subject_info.get("sample_processing_id", None),
                                             sample_id=subject_info.get("sample_id"),
                                             sample_type=subject_info.get("sample_type"),
                                             tissue=Ontology(id=subject_info.get("tissue_id", ""), lable=subject_info.get("tissue_label", "")),
                                             anatomic_site=subject_info.get("anatomic_site"),
                                             disease_state_sample=subject_info.get("disease_state_sample"),
                                             collection_time_point_relative=subject_info.get("collection_time_point_relative"),
                                             collection_time_point_relative_unit=Ontology(id=subject_info.get("collection_time_point_relative_unit_id", ""), lable=subject_info.get("collection_time_point_relative_unit_label", "")),
                                             collection_time_point_reference=subject_info.get("collection_time_point_reference"),
                                             biomaterial_provider=subject_info.get("biomaterial_provider"),
                                             tissue_processing=subject_info.get("tissue_processing"),
                                             cell_subset=Ontology(id=subject_info.get("cell_subset_id", ""), lable=subject_info.get("cell_subset_label", "")),
                                             cell_phenotype=subject_info.get("cell_phenotype"),
                                             cell_species=Ontology(id=subject_info.get("cell_species_id", ""), lable=subject_info.get("cell_species_label", "")),
                                             single_cell=str_to_bool((subject_info.get("single_cell"))),
                                             cell_number=subject_info.get("cell_number") if subject_info.get("cell_number") != '' else 0,
                                             cells_per_reaction=subject_info.get("cells_per_reaction") if subject_info.get("cells_per_reaction") != '' else 0,
                                             cell_storage=str_to_bool(subject_info.get("cell_storage")),
                                             cell_quality=subject_info.get("cell_quality"),
                                             cell_isolation=subject_info.get("cell_isolation"),
                                             cell_processing_protocol=subject_info.get("cell_processing_protocol"),
                                             template_class=TemplateClass(subject_info.get("template_class").upper()),
                                             template_quality=subject_info.get("template_quality"),
                                             template_amount=subject_info.get("template_amount"),
                                             template_amount_unit=Ontology(id=subject_info.get("template_amount_unit_id", ""), lable=subject_info.get("template_amount_unit_label", "")),
                                             library_generation_method=library_generation_method,
                                             library_generation_protocol=subject_info.get("library_generation_protocol"),
                                             library_generation_kit_version=subject_info.get("library_generation_kit_version"),
                                             pcr_target=subject_info.get("pcr_target"),
                                             complete_sequences=create_complete_sequences_enum(subject_info),
                                             physical_linkage=PhysicalLinkage(subject_info.get("physical_linkage")),
                                             sequencing_run_id=subject_info.get("sequencing_run_id"),
                                             total_reads_passing_qc_filter=subject_info.get("total_reads_passing_qc_filter"),
                                             sequencing_platform=subject_info.get("sequencing_platform"),
                                             sequencing_facility=subject_info.get("sequencing_facility"),
                                             sequencing_run_date=subject_info.get("sequencing_run_date"),
                                             sequencing_kit=subject_info.get("sequencing_kit"),
                                             sequencing_files=create_sequencing_data_object(subject_info)
                                             )
    
    sample_processing_list.append(sample_processing_obj)

    return sample_processing_list

def create_complete_sequences_enum(subject_info):
    """Create a CompleteSequences enum from subject information."""
    try:
        return CompleteSequences(subject_info.get("complete_sequences"))
    except:
        return CompleteSequences("partial")


def create_sequencing_data_object(subject_info):
    """Create a SequencingData object from subject information."""
    subject_info = fill_missing_required_fields(SequencingData,  subject_info)

    sequencing_data_obj = SequencingData(sequencing_data_id=subject_info.get("sequencing_data_id") if subject_info.get("sequencing_data_id") is not None else "",
                                         file_type=FileType(__root__=None),
                                         filename = subject_info.get("filename"),
                                         read_direction = ReadDirection(__root__=None),
                                         read_length = subject_info.get("read_length"),
                                         paired_filename = subject_info.get("paired_filename"),
                                         paired_read_direction = PairedReadDirection(__root__=None),
                                         paired_read_length = subject_info.get("paired_read_length"),
                                         index_filename = subject_info.get("index_filename"),
                                         index_length = subject_info.get("index_length"),
                                        )
    
    return sequencing_data_obj


def str_to_bool(value):
    """Convert a string to a boolean."""
    if isinstance(value, str):
        return value.lower() == 'true'
    return bool(value)


def create_subject_objects(subject_info):
    """Create Subject objects from subject information."""
    subject_info = fill_missing_required_fields(Subject,  subject_info)

    subject_object = Subject(subject_id=subject_info.get("subject_id"),
                             synthetic=subject_info.get("synthetic"),
                             species=Ontology(id=subject_info.get("species_id", ""), lable=subject_info.get("species_label", "")),
                             organism=Ontology(id=subject_info.get("organism_id", ""), lable=subject_info.get("organism_label", "")), 
                             sex=Sex(__root__=create_Sex_Enum(subject_info)),
                             age_min=subject_info.get("age_min"),
                             age_max=subject_info.get("age_max"),
                             age_unit=Ontology(id=subject_info.get("age_unit_id", ""), lable=subject_info.get("age_unit_label", "")),
                             age_event=subject_info.get("age_event"),
                             age=subject_info.get("age"),
                             ancestry_population=subject_info.get("ancestry_population"),
                             ethnicity=subject_info.get("ethnicity"),
                             race=subject_info.get("race", None),
                             strain_name=subject_info.get("strain_name", None),
                             linked_subjects=subject_info.get("linked_subjects", None),
                             link_type=subject_info.get("link_type", None),
                             diagnosis=create_diagnosis_list(subject_info),
                             genotype=create_subject_genotype_obj(subject_info))

    return subject_object


def create_Sex_Enum(subject_info):
    """Create a Sex enum from subject information."""
    try:
        return SexEnum(subject_info.get("sex"))
    except:
        return None


def create_subject_genotype_obj(subject_info):
    """Create a SubjectGenotype object from subject information."""
    subject_genotype_obj = SubjectGenotype(receptor_genotype_set=GenotypeSet(receptor_genotype_set_id=subject_info.get("receptor_genotype_set_id") if subject_info.get("receptor_genotype_set_id") is not None else "",
                                                                             genotype_class_list=create_genotype_model_list(subject_info)))
    return subject_genotype_obj


def create_genotype_model_list(subject_info):
    """Create a list of GenotypeModel objects from subject information."""
    # genotype_model_list = []
    # genotype_model_obj = GenotypeModel(receptor_genotype_id=)
    return None


def create_diagnosis_list(subject_info):
    """Create a list of Diagnosis objects from subject information."""
    subject_info = fill_missing_required_fields(Diagnosis,  subject_info)

    diagnosis_list = []
    diagnosis_object = Diagnosis(study_group_description=subject_info.get("study_group_description"),
                                 disease_diagnosis=Ontology(id=subject_info.get("disease_diagnosis_id", ""), lable=subject_info.get("disease_diagnosis_label", "")),
                                 disease_length=subject_info.get("disease_length"),
                                 disease_stage=subject_info.get("disease_stage"),
                                 prior_therapies=subject_info.get("prior_therapies"),
                                 immunogen=subject_info.get("immunogen"),
                                 intervention=subject_info.get("intervention"),
                                 medical_history=subject_info.get("medical_history"))
    diagnosis_list.append(diagnosis_object)
    return diagnosis_list


def create_study_object(subject_info):
    """Create a Study object from subject information."""
    subject_info = fill_missing_required_fields(Study,  subject_info)

    study_object = Study(study_id=subject_info.get("study_id"),
                         study_title=subject_info.get("study_title"),
                         study_type=Ontology(id=subject_info.get("study_type_id"), label=subject_info.get("study_type_lable")),
                         study_description=subject_info.get("study_description", None),
                         inclusion_exclusion_criteria=subject_info.get("inclusion_exclusion_criteria", ""),
                         grants=subject_info.get("grants", ""),
                         study_contact=subject_info.get("study_contact", None),
                         collected_by=subject_info.get("collected_by", ""),
                         lab_name=subject_info.get("lab_name", ""),
                         lab_address=subject_info.get("lab_address", ""),
                         submitted_by=subject_info.get("submitted_by", ""),
                         pub_ids=subject_info.get("pub_ids", ""),
                         keywords_study=create_keyword_study_list(subject_info), ###
                         adc_publish_date=create_date(subject_info, "adc_publish_date"),
                         adc_update_date=create_date(subject_info, "adc_update_date"),
                         )
    return study_object


def create_keyword_study_list(subject_info):
    """Create a list of KeywordsStudyEnum from subject information."""
    return [KeywordsStudyEnum('contains_ig')]


def create_date(subject_info, field):
    """Create a datetime object from subject information and field."""
    if subject_info.get(field) == "":
        return None
    
    else:
        return datetime.now()


def get_genomic_subjet_info(species, dataset, sample_id):
    """ Returns information on the selected sample """

    db = get_genomic_db(species, dataset)

    if db is None:
        return 'Bad species or dataset name', False

    sample = db.session.query(SAMPLE)\
        .filter(SAMPLE.sample_name == sample_id)\
        .one_or_none()

    if sample is None:
        return 'Bad sample name', False

    attribute_query = []

    for col in genomic_sample_filters.keys():
        if genomic_sample_filters[col]['field'] is not None:
            attribute_query.append(genomic_sample_filters[col]['field'])

    info = db.session.query(*attribute_query)\
        .filter(SAMPLE.sample_name == sample_id)\
        .join(Patient, Patient.id == SAMPLE.patient_id)\
        .join(SeqProtocol, SeqProtocol.id == SAMPLE.seq_protocol_id)\
        .join(TissuePro, TissuePro.id == SAMPLE.tissue_pro_id)\
        .join(DataPro, DataPro.id == SAMPLE.data_pro_id) \
        .join(STD, SAMPLE.study_id == STD.id)\
        .one_or_none()

    if info is not None:
        info = info._asdict()
        for k, v in info.items():
            if isinstance(v, datetime):
                info[k] = v.date().isoformat()

    return info, True


def get_airr_subjet_info(species, dataset, sample):
    """ Returns information on the selected sample """

    if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
        return 'Bad species or dataset name', False

    session = vdjbase_dbs[species][dataset].session
    attribute_query = []

    for col in sample_info_filters.keys():
        if sample_info_filters[col]['field'] is not None:
            attribute_query.append(sample_info_filters[col]['field'])

    info = session.query(*attribute_query)\
        .join(Airr_GenoDetection, Airr_GenoDetection.id == Airr_Sample.geno_detection_id)\
        .join(Airr_Patient, Airr_Patient.id == Airr_Sample.patient_id)\
        .join(Airr_SeqProtocol, Airr_SeqProtocol.id == Airr_Sample.seq_protocol_id)\
        .join(Airr_TissuePro, Airr_TissuePro.id == Airr_Sample.tissue_pro_id)\
        .join(Airr_DataPro, Airr_DataPro.id == Airr_Sample.data_pro_id) \
        .join(Airr_Study, Airr_Sample.study_id == Airr_Study.id)\
        .filter(Airr_Sample.sample_name == sample).one_or_none()

    if info:
        info = info._asdict()

        for k,v in info.items():
            if v:
                if isinstance(v, (dt.datetime, dt.date)):
                    info[k] = v.isoformat()

        haplotypes = session.query(HaplotypesFile.by_gene_s).join(SamplesHaplotype).join(Airr_Sample).filter(Airr_Sample.sample_name==sample).order_by(HaplotypesFile.by_gene_s).all()
        info['haplotypes'] = [(h[0]) for h in haplotypes]

    return info, True


def get_default_value(field_type: Any) -> Any:
    """
    Get the default value for a given field type.
    
    Args:
        field_type: The type of the field.
    
    Returns:
        The default value for the field type.
    """
    if get_origin(field_type) is Union:
        args = get_args(field_type)
        field_type = args[0] if args[1] is type(None) else args[1]

    if field_type == 'int':
        return 0
    elif field_type == 'float':
        return 0.0
    elif field_type == 'str':
        return ""
    elif field_type == 'bool':
        return False
    elif field_type == 'list':
        return []
    elif field_type == 'dict':
        return {}
    elif field_type == 'datetime':
        return datetime.now()
    elif field_type == 'date':
        return datetime.now()
    else:
        if issubclass(globals()[field_type] , Enum):
            # Get the first value of the Enum
            return next(iter(globals()[field_type])).value

    return None


def fill_missing_required_fields(model_cls: BaseModel, data: dict) -> dict:
    """
    Fill missing required fields in the given data with default values.
    
    Args:
        model_cls: The Pydantic model class.
        data: The data dictionary.
    
    Returns:
        The data dictionary with missing required fields filled.
    """
    filled_data = data.copy()

    for field_name, field_info in model_cls.__fields__.items():
        if field_name in data:
            print(field_name)
            if field_name == 'sequencing_run_date':
                breakpoint()
            if is_required(field_info):
                field_type = model_cls.__annotations__[field_name]
                if data[field_name] is None or (field_type != 'str' and data[field_name] == ''):
                    if field_type != 'List[KeywordsStudyEnum]':
                        default_value = get_default_value(field_type)
                    if default_value is not None:
                        filled_data[field_name] = default_value

    return filled_data

def is_required(field_info: FieldInfo) -> bool:
    """
    Check if a field is required.
    
    Args:
        field_info: The field information.
    
    Returns:
        Boolean indicating if the field is required.
    """
    if 'required=True' in str(field_info):
        return True
    
    return False
