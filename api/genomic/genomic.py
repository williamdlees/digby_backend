# Services related to genomic sequences and features
import os

from flask import request
from flask_restx import Resource, reqparse

from api.reports.genotypes import process_genomic_genotype
from api.restx import api
from sqlalchemy import inspect, func, distinct
from math import ceil
from werkzeug.exceptions import BadRequest

from api.system.system import digby_protected
from db.genomic_db import RefSeq, Feature, Sequence, SampleSequence, Gene
from db.genomic_airr_model import Sample, Study, Subject, SeqProtocol, TissuePro, DataPro, Base
from db.genomic_api_query_filters import genomic_sequence_filters, genomic_subject_filters

import json
from datetime import datetime
from sqlalchemy_filters import apply_filters
from decimal import Decimal
from app import app, genomic_dbs
from api.vdjbase.vdjbase import get_vdjbase_species


# Return SqlAlchemy row as a dict, using correct column names
def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


GENOMIC_SAMPLE_PATH = os.path.join(app.config['STATIC_PATH'], 'study_data/Genomic/samples')

ns = api.namespace('genomic', description='Genomic data and annotations')


def get_genomic_species():
    return list(genomic_dbs.keys())


@ns.route('/species')
@api.response(404, 'No species available!')
class SpeciesApi(Resource):
    @digby_protected()
    def get(self):
        """ Returns the list of species for which information is held """
        genomic_sp = get_genomic_species()
        vdjbase_sp = get_vdjbase_species()
        return list(set(genomic_sp) | set(vdjbase_sp))


def get_genomic_datasets(species):
    sp = get_genomic_species()

    if species in sp:
        datasets = [{'dataset': d, 'locus': d} for d in genomic_dbs[species].keys() if 'description' not in d]
        return datasets
    else:
        return []


@ns.route('/data_sets/<string:species>')
@api.response(404, 'No data sets available for that species.')
class DataSetAPI(Resource):
    @digby_protected()
    def get(self, species):
        """ Returns the list of data sets for the selected species """
        return get_genomic_datasets(species)


def get_genomic_db(species, dataset):
    if species in genomic_dbs and dataset in genomic_dbs[species]:
        return genomic_dbs[species][dataset]
    else:
        return None


@ns.route('/subject_info/<string:species>/<string:dataset>/<string:sample_id>')
class SubjectInfoApi(Resource):
    @digby_protected()
    def get(self, species, dataset, sample_id):
        """ Returns information on the selected sample """

        db = get_genomic_db(species, dataset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        sample = db.session.query(Sample)\
            .filter(Sample.identifier == sample_id)\
            .one_or_none()

        if sample is None:
            raise BadRequest('Bad sample name')

        attribute_query = []

        for col in genomic_subject_filters.keys():
            if genomic_subject_filters[col]['field'] is not None:
                attribute_query.append(genomic_subject_filters[col]['field'])

        info = db.session.query(*attribute_query)\
            .join(Subject, Subject.id == Sample.subject_id)\
            .join(Study, Study.id == Subject.study_id)\
            .filter(Sample.identifier == sample_id)\
            .one_or_none()

        if info is not None:
            info = info._asdict()
            for k, v in info.items():
                if isinstance(v, datetime):
                    info[k] = v.date().isoformat()

        return info


range_arguments = reqparse.RequestParser()
range_arguments.add_argument('start', type=int, required=True, location='args')
range_arguments.add_argument('end', type=int, required=True, location='args')


def enumerate_feature(f):
    ret = {
        'type': f.feature,
#        'name' : f.name,
        'start': f.start,
        'end': f.end,
        'score': f.score if f.score else 0,
        'strand': 1 if f.strand == '+' else 0,
        'phase': f.frame if f.frame else 0,
    }

    if f.feature_id:
        ret['uniqueID'] = str(f.feature_id)

    if f.feature == 'gene':
        ret['name'] = f.name

    return ret

'''
genomic_sequence_filters = {
    'name': {'model': 'Sequence', 'field': Sequence.name, 'sort': 'gene'},
    'imgt_name': {'model': 'Sequence', 'field': Sequence.imgt_name},
    'type': {'model': 'Sequence', 'field': Sequence.type},
    'novel': {'model': 'Sequence', 'field': Sequence.novel},
    'deleted': {'model': 'Sequence', 'field': Sequence.deleted},
    'functional': {'model': 'Sequence', 'field': Sequence.functional},
    'notes': {'model': 'Sequence', 'field': Sequence.notes},
    'sequence': {'model': 'Sequence', 'field': Sequence.sequence, 'no_uniques': True},
    'gapped_sequence': {'model': 'Sequence', 'field': Sequence.gapped_sequence, 'no_uniques': True},
    'appearances': {'model': 'Sequence', 'field': Sequence.appearances, 'fieldname': 'appearances', 'sort': 'numeric'},

    'gene_name': {'model': 'Gene', 'field': Gene.name.label('gene_name'), 'sort': 'gene', 'fieldname': 'gene_name`'},

    'sample_identifier': {'model': None, 'fieldname': 'sample_identifier'},
    'dataset': {'model': None, 'fieldname': 'dataset'},
}
'''

genomic_sequence_bool_values = {
    'novel': ('Novel', '(blank)'),
    'deleted': ('Deleted', '(blank)'),
}

filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('page_number', type=int, location='args')
filter_arguments.add_argument('page_size', type=int, location='args')
filter_arguments.add_argument('filter', type=str, location='args')
filter_arguments.add_argument('sort_by', type=str, location='args')
filter_arguments.add_argument('cols', type=str, location='args')

# borrowed from sqlalchemy-filters


OPERATORS = {
    '==': lambda f, a: f == a,
    '!=': lambda f, a: f != a,
    '>': lambda f, a: f > a,
    '<': lambda f, a: f < a,
    '>=': lambda f, a: f >= a,
    '<=': lambda f, a: f <= a,
    'like': lambda f, a: f.like(a),
    'ilike': lambda f, a: f.ilike(a),
    'not_ilike': lambda f, a: ~f.ilike(a),
    'in': lambda f, a: f.in_(a),
    'not_in': lambda f, a: ~f.in_(a),
    'any': lambda f, a: f.any(a),
    'not_any': lambda f, a: func.not_(f.any(a)),
}


@ns.route('/sequences/<string:species>/<string:genomic_datasets>')
@api.response(404, 'Reference sequence not found.')
class SequencesAPI(Resource):
    @digby_protected()
    @api.expect(filter_arguments, validate=True)
    def get(self, species, genomic_datasets):
        """ Returns nucleotide sequences from selected reference or multiple references (separate multiple reference names with ',')  """
        args = filter_arguments.parse_args(request)

        required_cols = json.loads(args['cols'])
        genomic_datasets = genomic_datasets.split(',')
        ret = find_genomic_sequences(required_cols, genomic_datasets, species, json.loads(args['filter']) if args['filter'] else [])

        gene_order = {}
        set_index = 0
        for ds in genomic_datasets:
            db = get_genomic_db(species, ds)
            go = db.session.query(Gene.name, Gene.alpha_order).all()
            for gene, order in go:
                gene_order[gene] = int(order) + set_index

        uniques = {}
        for f in required_cols:
            uniques[f] = []
        uniques['dataset'] = genomic_datasets

        for s in ret:
            for f in required_cols:
                if 'field' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['field'] is not None and 'no_uniques' not in genomic_sequence_filters[f]:
                    el = s[f]
                    if isinstance(el, datetime):
                        el = el.date().isoformat()
                    elif isinstance(el, Decimal):
                        el = int(el)
                    elif isinstance(el, str) and len(el) == 0:
                        el = '(blank)'
                    elif isinstance(el, bool):
                        if f in genomic_sequence_bool_values:
                            el = genomic_sequence_bool_values[f][0 if el else 1]
                    if el not in uniques[f]:
                        uniques[f].append(el)


        def allele_sort_key(name):
            if gene_order is None or len(gene_order) == 0:
                return fallback_allele_sort_key(name)

            if '*' in name:
                gene = name.split('*')
            else:
                gene = (name, '')

            return((gene_order[gene[0]] if gene[0] in gene_order else 999, gene[1]))


        def fallback_allele_sort_key(x):
            if x is None or x == '' or 'IG' not in x:
                return ''
            name = x.split('IGH')[1]
            seg = name[:1]
            name = name[1:]
            allele = '0000'
            if '*' in name:
                name, allele = name.split('*')
                allele = allele.replace('0', '').zfill(4)
            if '.' in name:
                fam, num = name.split('.', 1)
            elif '-' in name:
                fam, num, = name.split('-', 1)
            else:
                fam = ''
                num = name

            if '.' in num:
                num = num.split('.', 1)[0]
            if '-' in num:
                num = num.split('-', 1)[0]
            num = num.zfill(8)
            return seg+num+allele


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
                if 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'numeric':
                    uniques[f].sort(key=num_sort_key)
                elif 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'gene':
                    uniques[f].sort(key=allele_sort_key)
                else:
                    uniques[f].sort(key=lambda x: (x is None or x == '', x))
            except:
                pass

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None) else []
        if len(sort_specs) == 0:
            sort_specs = [{'field': 'name', 'order': 'asc'}]

        for spec in sort_specs:
            if spec['field'] in genomic_sequence_filters.keys():
                f = spec['field']
                if 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'numeric':
                    ret = sorted(ret, key=lambda x: num_sort_key(x[spec['field']]), reverse=(spec['order'] == 'desc'))
                elif 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'gene':
                    ret = sorted(ret, key=lambda x: allele_sort_key(x[spec['field']]), reverse=(spec['order'] == 'desc'))
                else:
                    ret = sorted(ret, key=lambda x : (x[spec['field']] is None,  x[spec['field']]), reverse=(spec['order'] == 'desc'))

        total_size = len(ret)

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first : first + args['page_size']]

        return {
            'sequences': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
        }


def find_genomic_sequences(required_cols, genomic_datasets, species, genomic_filters):
    ret = []
    for dataset in genomic_datasets:
        db = get_genomic_db(species, dataset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        for col in required_cols:
            if col not in genomic_sequence_filters.keys():
                raise BadRequest('Bad column string %s' % col)

        if 'sequence' in required_cols and 'gapped_sequence' not in required_cols:
            required_cols.append('gapped_sequence')

        attribute_query = [
            genomic_sequence_filters['name']['field']]  # the query requires the first field to be from Sequence

        for col in required_cols:
            if col != 'name' and 'field' in genomic_sequence_filters[col] and genomic_sequence_filters[col]['field'] is not None:
                attribute_query.append(genomic_sequence_filters[col]['field'])

        seq_query = db.session.query(*attribute_query)

        if 'gene_name' in required_cols:
            seq_query = seq_query.join(Gene, Sequence.gene_id == Gene.id)

        filter_spec = []
        sample_count_filters = []
        sample_id_filter = None

        if len(genomic_filters) > 0:
            for f in genomic_filters:
                try:
                    if 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'sample_count':
                        sample_count_filters.append(f)
                    elif 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'sample_identifier':
                        sample_id_filter = f
                    elif f['field'] == 'dataset':
                        if f['op'] == 'in' and dataset not in f['value']:
                            continue  # just going to ignore other criteria I'm afraid
                    else:
                        f['model'] = genomic_sequence_filters[f['field']]['model']
                        if 'fieldname' in genomic_sequence_filters[f['field']]:
                            f['field'] = genomic_sequence_filters[f['field']]['fieldname']
                        if f['field'] in genomic_sequence_bool_values:
                            value = []
                            for v in f['value']:
                                value.append('1' if v == genomic_sequence_bool_values[f['field']][0] else '0')
                            f['value'] = value
                        elif '(blank)' in f['value']:
                            value_specs = [
                                {'model': genomic_sequence_filters[f['field']]['model'], 'field': f['field'], 'op': 'is_null', 'value': ''},
                                {'model': genomic_sequence_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': ''},
                            ]

                            for v in f['value']:
                                if v != '(blank)':
                                    value_specs.append({'model': genomic_sequence_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': v})
                            
                            f = {'or': value_specs}

                        filter_spec.append(f)
                except Exception as e:
                    raise BadRequest(f'Bad filter string: {f}: {e}')

        seq_query = apply_filters(seq_query, filter_spec)

        for f in sample_count_filters:
            if f['op'] in OPERATORS:
                seq_query = seq_query.having(OPERATORS[f['op']](func.count(Subject.name), f['value']))

        appears = {}

        if sample_id_filter is not None:
            filtered_sample_ids = []

            names_to_ids = {}
            for name, id in db.session.query(Sample.identifier, Sample.id).all():
                names_to_ids[name] = id

            for names in sample_id_filter['value'].items():
                if names[0] == dataset:
                    for n in names[1]:
                        filtered_sample_ids.append(names_to_ids[n])

            if not filtered_sample_ids:
                continue

            ids_to_names = {}
            for id, name in  db.session.query(Sequence.id, Sequence.name).all():
                ids_to_names[id] = name

            sequence_id_query = db.session.query(Sequence.id.distinct()) \
                .join(SampleSequence, Sequence.id == SampleSequence.sequence_id) \
                .join(Sample, Sample.id == SampleSequence.sample_id)

            filtered_sequence_ids = sequence_id_query.filter(Sample.id.in_(filtered_sample_ids)).all()
            filtered_sequence_ids = [x[0] for x in filtered_sequence_ids]
            seq_query = seq_query.filter(Sequence.id.in_(filtered_sequence_ids))

            subseqs = db.session.query(SampleSequence.sequence_id, SampleSequence.sample_id)\
                .filter(SampleSequence.sample_id.in_(filtered_sample_ids))\
                .filter(SampleSequence.sequence_id.in_(filtered_sequence_ids))\
                .all()

            appearances = {}

            for sequence_id, subject_id in subseqs:
                if sequence_id not in appearances:
                    appearances[sequence_id] = []
                appearances[sequence_id].append(subject_id)

            for k, v in appearances.items():
                appears[ids_to_names[k]] = len(set(v))

        seqs = seq_query.all()

        for r in seqs:
            s = r._asdict()

            if len(appears):
                if s['name'] in appears:
                    s['appearances'] = appears[s['name']]
                else:
                    s['appearances'] = 0

            for k, v in s.items():
                if isinstance(v, datetime):
                    s[k] = v.date().isoformat()
                elif isinstance(v, Decimal):
                    s[k] = int(v)
            s['dataset'] = dataset

            ret.append(s)

    return ret


@ns.route('/feature_pos/<string:species>/<string:dataset>/<string:ref_seq_name>/<string:feature_string>')
@api.response(404, 'Reference sequence not found.')
class FeaturePosAPI(Resource):
    @digby_protected()
    def get(self, species, dataset, ref_seq_name, feature_string):
        """ Returns the position of the first feature matching the specified string """

        db = get_genomic_db(species, dataset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        ref_seq = db.session.query(RefSeq).filter(RefSeq.name == ref_seq_name).one_or_none()

        if ref_seq is None:
            raise BadRequest('Bad ref_seq name')
        
        if feature_string.upper() == ref_seq_name.upper():
            return [{'chromosome': ref_seq_name, 'start': 1, 'end': ref_seq.length}]

        features = db.session.query(Feature).join(RefSeq).filter(RefSeq.id == ref_seq.id).filter(Feature.name.contains(feature_string)).all()

        if features:
            start = 99999999
            end = 0

            for feature in features:
                start = min(start, feature.start)
                end = max(end, feature.end)

            return[{'chromosome': ref_seq_name, 'start': start, 'end': end}]
        else:
            return []

'''
genomic_subject_filters = {
    'subject_identifier': {'model': 'Subject', 'field': Subject.identifier.label('subject_identifier'), 'fieldname': 'identifier', 'sort': 'underscore'},
    'name_in_study': {'model': 'Subject', 'field': Subject.name_in_study},
    'mother_in_study': {'model': 'Subject', 'field': Subject.mother_in_study},
    'father_in_study': {'model': 'Subject', 'field': Subject.father_in_study},
    'age': {'model': 'Subject', 'field': Subject.age},
    'sex': {'model': 'Subject', 'field': Subject.sex},
    'self_ethnicity': {'model': 'Subject', 'field': Subject.self_ethnicity},
    'grouped_ethnicity': {'model': 'Subject', 'field': Subject.grouped_ethnicity},
    'population': {'model': 'Subject', 'field': Subject.population},
    'population_abbr': {'model': 'Subject', 'field': Subject.population_abbr},
    'super_population': {'model': 'Subject', 'field': Subject.super_population},
  
    'id': {'model': 'Sample', 'field': Sample.id},
    'sample_name_in_study': {'model': 'Sample', 'field': Sample.name_in_study.label('sample_name_in_study'), 'fieldname': 'name_in_study'},
    'sample_identifier': {'model': 'Sample', 'field': Sample.identifier.label('sample_identifier'), 'fieldname': 'identifier', 'sort': 'underscore'},
    'annotation_path': {'model': 'Sample', 'field': Sample.annotation_path},
    'annotation_method': {'model': 'Sample', 'field': Sample.annotation_method},
    'annotation_format': {'model': 'Sample', 'field': Sample.annotation_format},
    'annotation_reference': {'model': 'Sample', 'field': Sample.annotation_reference},
    'locus_coverage': {'model': 'Sample', 'field': Sample.locus_coverage},
    'sequencing_platform': {'model': 'Sample', 'field': Sample.sequencing_platform},
    'assembly_method': {'model': 'Sample', 'field': Sample.assembly_method},
    'DNA_source': {'model': 'Sample', 'field': Sample.DNA_source},

    'study_name': {'model': 'Study', 'field': Study.study_name.label('study_name'), 'fieldname': 'study_name'},
    'study_title': {'model': 'Study', 'field': Study.title.label('study_title'), 'fieldname': 'study_title'},
    'study_id': {'model': 'Study', 'field': Study.study_id.label('study_id'), 'fieldname': 'study_id'},
    'study_date': {'model': 'Study', 'field': Study.date.label('study_date'), 'fieldname': 'date'},
    'study_description': {'model': 'Study', 'field': Study.description.label('study_description'), 'fieldname': 'description'},
    'institute': {'model': 'Study', 'field': Study.institute},
    'researcher': {'model': 'Study', 'field': Study.researcher},
    'reference': {'model': 'Study', 'field': Study.reference},
    'contact': {'model': 'Study', 'field': Study.contact},

    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset'},
}
'''

@ns.route('/subjects/<string:species>/<string:genomic_datasets>')
@api.response(404, 'Reference sequence not found.')
class SubjectsAPI(Resource):
    @digby_protected()
    @api.expect(filter_arguments, validate=True)
    def get(self, species, genomic_datasets):
        """ Returns a list of subjects in the selected datasets  """
        args = filter_arguments.parse_args(request)

        required_cols = json.loads(args['cols']) if 'cols' in args and args['cols'] else list(genomic_subject_filters.keys())

        for col in required_cols:
            if col not in genomic_subject_filters.keys():
                raise BadRequest('Bad filter string %s' % args['filter'])

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
            ret = ret[first : first + args['page_size']]

        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
        }


def find_genomic_samples(attribute_query, species, genomic_datasets, genomic_filters):
    results = []
    for dataset in genomic_datasets:
        db = get_genomic_db(species, dataset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        sample_query = db.session.query(*attribute_query)\
            .select_from(Sample)\
            .join(Subject, Sample.subject_id == Subject.id)\
            .join(Study, Study.id == Subject.study_id)

        allele_filters = None

        filter_spec = []
        for f in genomic_filters:
            try:
                if f['field'] == 'allele':
                    allele_filters = f
                elif f['field'] == 'dataset':
                    if f['op'] == 'in' and dataset not in f['value']:
                        continue        # just going to ignore other criteria I'm afraid
                else:
                    f['model'] = genomic_subject_filters[f['field']]['model']
                    if 'fieldname' in genomic_subject_filters[f['field']]:
                        f['field'] = genomic_subject_filters[f['field']]['fieldname']
                    if '(blank)' in f['value']:
                        value_specs = [
                            {'model': genomic_subject_filters[f['field']]['model'], 'field': f['field'], 'op': 'is_null', 'value': ''},
                            {'model': genomic_subject_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': ''},
                        ]

                        for v in f['value']:
                            if v != '(blank)':
                                value_specs.append({'model': genomic_subject_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': v})
                        
                        f = {'or': value_specs}

                    filter_spec.append(f)
            except Exception as e:
                raise BadRequest(f'Bad filter string: {f}: {e}')

        if len(filter_spec) > 0:
            sample_query = apply_filters(sample_query, filter_spec)

        if allele_filters is not None:
            samples_with_alleles = db.session.query(Sample.identifier)\
                .join(Subject, Sample.subject_id == Subject.id)\
                .join(SampleSequence, SampleSequence.sample_id == Sample.id)\
                .join(Sequence, SampleSequence.sequence_id == Sequence.id)\
                .filter(Sequence.name.in_(allele_filters['value']))\
                .all()
            samples_with_alleles = [x[0] for x in samples_with_alleles]
            sample_query = sample_query.filter(Sample.identifier.in_(samples_with_alleles))

        samples = sample_query.all()

        for s in samples:
            r = s._asdict()
            for k in list(r.keys()):
                v = r[k]
                if isinstance(v, datetime):
                    r[k] = v.date().isoformat()
                elif k == 'annotation_path':
                    if v is None:
                        app.logger.error('No annotation path for sample %s' % r['sample_identifier'])
                        r[k] = ''
                    elif 'http' not in v:
                        r[k] = os.path.join(app.config['STATIC_LINK'], 'study_data/Genomic/samples', r[k])
            r['dataset'] = dataset
            results.append(r)

    return results


def find_genomic_filter_params(species, genomic_datasets):
    genes = []
    gene_types = []

    for dset in genomic_datasets:
        db = get_genomic_db(species, dset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        g = db.session.query(Gene.name).all()
        genes.extend(g)
        g_t = db.session.query(Gene.type).distinct().all()
        gene_types.extend(g_t)

    genes = sorted(set([gene[0] for gene in genes]))
    gene_types = sorted(set([gene_type[0] for gene_type in gene_types]))

    params = [
        {
          "id": "f_pseudo_genes",
          "type": "boolean",
          "label": "Include pseudogenes"
        },
    ]

    params.extend([
        {
            "id": "f_gene_types",
            "type": "multi_select",
            "label": "Only process selected gene types (leave blank for all)",
            "options": gene_types
        },
        {
            "id": "f_genes",
            "type": "multi_select",
            "label": "Only process selected genes (leave blank for all)",
            "options": genes

        }
    ])

    return params


@ns.route('/assemblies/<string:species>/<string:data_sets>')
class AssemblyAPI(Resource):
    @digby_protected()
    def get(self, species, data_sets):
        """ Returns the list of annotated assemblies for the selected species and datasets """
        results = []
        for dataset in data_sets.split(','):
            db = get_genomic_db(species, dataset)

            if db is None:
                raise BadRequest('Bad species or dataset name')

            assemblies = db.session.query(distinct(RefSeq.name)).all()
            results.extend(assemblies)

        return [{'assembly': row[0]} for row in results]


@ns.route('/genotype/<string:species>/<string:subject_name>')
class SamplesApi(Resource):
    @digby_protected()
    def get(self, species, subject_name):
        """ Returns the inferred genotype (in MiAIRR format) of the specified subject """

        if species not in genomic_dbs:
            return None, 404

        genotypes = []

        for dataset in genomic_dbs[species].keys():
            if '_description' not in dataset:
                genotype = self.single_genotype(species, dataset, subject_name)
                if genotype:
                    genotypes.append(genotype)

        if not genotypes:
            return None, 404

        ret = {
            'GenotypeSet': {
                'receptor_genotype_id': 'Genomic_genotype_set_' + subject_name,
                'genotype_class_list': genotypes
            }

        }

        return ret

    def single_genotype(self, species, dataset, subject_name):
        print(dataset)
        session = genomic_dbs[species][dataset].session
        sample = session.query(Subject).filter(Subject.identifier == subject_name).one_or_none()

        if not sample:
            return None

        reference_set_version = sample.reference_set_version
        genotype = process_genomic_genotype(subject_name, [], session, True, False)
        documented = []
        undocumented = []
        deleted = []
        for row in genotype.itertuples():
            gene_type = row.gene[3]

            for allele in row.GENOTYPED_ALLELES.split(','):
                allele_name = row.gene + '*' + allele
                res = session.query(Sequence.sequence, Sequence.novel).filter(Sequence.name == allele_name).one_or_none()
                if res:
                    seq, novel = res

                    if novel:
                        undocumented.append({'allele_name': allele_name, 'germline_set_ref': reference_set_version, 'sequence': seq, 'phasing': 0})
                    else:
                        documented.append({'allele_name': allele_name, 'germline_set_ref': reference_set_version, 'phasing': 0})
        ret = {
            'receptor_genotype_id': 'Tigger_genotype_' + subject_name + '_' + dataset,
            'locus': dataset,
            'documented_alleles': documented,
            'undocumented_alleles': undocumented,
            'deleted_genes': deleted,
            'inference_process': 'genomic_sequencing',
            'genotyping_tool': sample.annotation_method,
        }
        return ret
