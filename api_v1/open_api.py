import json
import datetime
from flask import Blueprint, request, jsonify, Response, send_from_directory
from api.vdjbase.vdjbase import get_vdjbase_species, find_datasets
from api.genomic.genomic import get_genomic_species, get_genomic_datasets, genomic_subject_filters, find_genomic_samples, ceil
find_datasets
from schema.models import Ontology, ErrorResponse, SpeciesResponse,Dataset, DatasetsResponse

api_bp = Blueprint('api_v1', __name__)

@api_bp.route('/<type>/species', methods=['GET'])
def get_species(type):
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
    if type == "genomic":
        subjects_list, is_currect = get_genomic_list_subjects(species, dataset)
        if not is_currect:
            error_response = ErrorResponse(message=str(subjects_list))
            return jsonify(error_response.dict()), 500
        
        else:
            pass

    elif type == "airrseq":
        data_sets = find_datasets(species)
    pass


def get_genomic_list_subjects(species, genomic_datasets):
    args = {
        "page_number": 0,
        "page_size": 25,
        "filter": None,
        "sort_by": None,
        "cols": '["subject_identifier","sample_identifier"]'
    }

    required_cols = json.loads(args['cols']) if 'cols' in args and args['cols'] else list(genomic_subject_filters.keys())

    for col in required_cols:
        if col not in genomic_subject_filters.keys():
            return 'Bad filter string %s' % args['filter'], False

    if 'study_name' not in required_cols:
        required_cols = ['study_name'] + required_cols
    if 'dataset' not in required_cols:
        required_cols.append('dataset')
    if 'sample_identifier' not in required_cols:
        required_cols.append('sample_identifier')
    if 'annotation_path' in required_cols:
        if 'annotation_method' not in required_cols:
            required_cols.append('annotation_method')
        if 'annotation_reference' not in required_cols:
            required_cols.append('annotation_reference')

    attribute_query = [genomic_subject_filters['id']['field']]        # the query requires the first field to be from Sample

    for col in required_cols:
        if col != 'id' and genomic_subject_filters[col]['field'] is not None:
            attribute_query.append(genomic_subject_filters[col]['field'])

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
            if genomic_subject_filters[f]['field'] is not None and 'no_uniques' not in genomic_subject_filters[f]:
                el = s[f]
                if isinstance(el, datetime):
                    el = el.date().isoformat()
                elif isinstance(el, str) and len(el) == 0:
                    el = '(blank)'
                if el not in uniques[f]:
                    uniques[f].append(el)
        if filter_applied:
            uniques['names_by_dataset'][s['dataset']].append(s['sample_identifier'])

    for f in required_cols:
        try:
            if 'sort' in genomic_subject_filters[f] and genomic_subject_filters[f]['sort'] == 'numeric':
                uniques[f].sort(key=num_sort_key)
            elif 'sort' in genomic_subject_filters[f] and genomic_subject_filters[f]['sort'] == 'underscore':
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
        if f in genomic_subject_filters.keys():
            if 'sort' in genomic_subject_filters[f] and genomic_subject_filters[f]['sort'] == 'underscore':
                ret = sorted(ret, key=lambda x: name_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
            elif 'sort' in genomic_subject_filters[f] and genomic_subject_filters[f]['sort'] == 'numeric':
                ret = sorted(ret, key=lambda x: num_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
            else:
                ret = sorted(ret, key=lambda x: ((x[f] is None or x[f] == ''),  x[f]), reverse=(spec['order'] == 'desc'))

    total_size = len(ret)

    if args['page_size']:
        first = (args['page_number']) * args['page_size']
        ret = ret[first: first + args['page_size']]

    return {
        'samples': ret,
        'uniques': uniques,
        'total_items': total_size,
        'page_size': args['page_size'],
        'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
    }
