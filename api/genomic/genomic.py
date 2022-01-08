# Services related to genomic sequences and features

from flask import request
from flask_restx import Resource, reqparse, inputs
from api.restx import api
from sqlalchemy import inspect, or_, func, distinct
from math import ceil
from werkzeug.exceptions import BadRequest

from api.system.system import digby_protected
from db.genomic_db import RefSeq, Feature, Sequence, Subject, Study, SubjectSequence, Gene
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

GENOMIC_SAMPLE_PATH = 'study_data/Genomic/samples'

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


@ns.route('/subject_info/<string:species>/<string:dataset>/<string:subject>')
class SubjectInfoApi(Resource):
    @digby_protected()
    def get(self, species, dataset, subject):
        """ Returns information on the selected sample """

        db = get_genomic_db(species, dataset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        sample = db.session.query(Subject)\
            .filter(Subject.identifier == subject)\
            .one_or_none()

        if sample is None:
            raise BadRequest('Bad subject name')

        attribute_query = []

        for col in genomic_subject_filters.keys():
            if genomic_subject_filters[col]['field'] is not None:
                attribute_query.append(genomic_subject_filters[col]['field'])

        info = db.session.query(*attribute_query)\
            .join(Study, Study.id == Subject.study_id)\
            .filter(Subject.identifier == subject)\
            .one_or_none()

        if info is not None:
            info = info._asdict()
            for k, v in info.items():
                if isinstance(v, datetime):
                    info[k] = v.date().isoformat()

        return info


range_arguments = reqparse.RequestParser()
range_arguments.add_argument('start', type=int, required=True)
range_arguments.add_argument('end', type=int, required=True)


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


genomic_sequence_filters = {
    'name': {'model': 'Sequence', 'field': Sequence.name, 'sort': 'gene'},
    'gene': {'model': 'Sequence', 'field': Sequence.gene, 'sort': 'gene'},
    'imgt_name': {'model': 'Sequence', 'field': Sequence.imgt_name},
    'type': {'model': 'Sequence', 'field': Sequence.type},
    'novel': {'model': 'Sequence', 'field': Sequence.novel},
    'deleted': {'model': 'Sequence', 'field': Sequence.deleted},
    'functional': {'model': 'Sequence', 'field': Sequence.functional},
    'sequence': {'model': 'Sequence', 'field': Sequence.sequence, 'no_uniques': True},
    'gapped_sequence': {'model': 'Sequence', 'field': Sequence.gapped_sequence, 'no_uniques': True},
    'appearances': {'model': 'Sequence', 'field': Sequence.appearances, 'fieldname': 'appearances', 'sort': 'numeric'},

    'subject_id': {'model': None, 'fieldname': 'subject_id'},
    'dataset': {'model': None, 'fieldname': 'dataset'},
}

genomic_sequence_bool_values = {
    'novel': ('Novel', '(blank)'),
    'deleted': ('Deleted', '(blank)'),
}

filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('page_number', type=int)
filter_arguments.add_argument('page_size', type=int)
filter_arguments.add_argument('filter', type=str)
filter_arguments.add_argument('sort_by', type=str)
filter_arguments.add_argument('cols', type=str)

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
        ret = find_genomic_sequences(required_cols, genomic_datasets, species, json.loads(args['filter']))

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

        for s in ret:
            for f in required_cols:
                if genomic_sequence_filters[f]['field'] is not None and 'no_uniques' not in genomic_sequence_filters[f]:
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
            'pages': ceil((total_size*1.0)/args['page_size'])
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
            if col != 'name' and genomic_sequence_filters[col]['field'] is not None:
                attribute_query.append(genomic_sequence_filters[col]['field'])

        seq_query = db.session.query(*attribute_query)

        filter_spec = []
        sample_count_filters = []
        subject_id_filter = None

        if len(genomic_filters) > 0:
            for f in genomic_filters:
                try:
                    if 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'sample_count':
                        sample_count_filters.append(f)
                    elif 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'subject_id':
                        subject_id_filter = f
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
                            f['value'].append('')
                        filter_spec.append(f)
                except:
                    raise BadRequest('Bad filter string %s' % f)

        seq_query = apply_filters(seq_query, filter_spec)

        for f in sample_count_filters:
            if f['op'] in OPERATORS:
                seq_query = seq_query.having(OPERATORS[f['op']](func.count(Subject.name), f['value']))

        if subject_id_filter is not None:
            filtered_subject_ids = []

            for names in subject_id_filter['value'].items():
                if names[0] == dataset:
                    subject_ids = db.session.query(Subject.id.distinct()) \
                        .filter(Subject.identifier.in_(names[1])).all()
                    subject_ids = [x[0] for x in subject_ids]
                    filtered_subject_ids.extend(subject_ids)

            if not filtered_subject_ids:
                continue

            sequence_id_query = db.session.query(Sequence.id.distinct()) \
                .join(SubjectSequence, Sequence.id == SubjectSequence.sequence_id) \
                .join(Subject, Subject.id == SubjectSequence.sample_id)

            filtered_sequence_ids = sequence_id_query.filter(Subject.id.in_(filtered_subject_ids)).all()
            filtered_sequence_ids = [x[0] for x in filtered_sequence_ids]
            seq_query = seq_query.filter(Sequence.id.in_(filtered_sequence_ids))

        seqs = seq_query.all()

        for r in seqs:
            s = r._asdict()
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

genomic_subject_filters = {
    'id': {'model': 'Subject', 'field': Subject.id},
    'identifier': {'model': 'Subject', 'field': Subject.identifier, 'sort': 'underscore'},
    'name_in_study': {'model': 'Subject', 'field': Subject.name_in_study},
    'age': {'model': 'Subject', 'field': Subject.age},
    'sex': {'sex': 'Subject', 'field': Subject.sex},
    'annotation_path': {'annotation_path': 'Subject', 'field': Subject.annotation_path},
    'annotation_method': {'annotation_method': 'Subject', 'field': Subject.annotation_method},
    'annotation_format': {'annotation_format': 'Subject', 'field': Subject.annotation_format},
    'annotation_reference': {'annotation_reference': 'Subject', 'field': Subject.annotation_reference},

    'study_name': {'model': 'Study', 'field': Study.name.label('study_name'), 'fieldname': 'study_name'},
    'study_date': {'model': 'Study', 'field': Study.date.label('study_date'), 'fieldname': 'study_date'},
    'study_description': {'model': 'Study', 'field': Study.description.label('study_description'), 'fieldname': 'study_description'},
    'institute': {'model': 'Study', 'field': Study.institute},
    'researcher': {'model': 'Study', 'field': Study.researcher},
    'reference': {'model': 'Study', 'field': Study.reference},
    'contact': {'model': 'Study', 'field': Study.contact},

    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset'},
}


@ns.route('/subjects/<string:species>/<string:genomic_datasets>')
@api.response(404, 'Reference sequence not found.')
class SubjectsAPI(Resource):
    @digby_protected()
    @api.expect(filter_arguments, validate=True)
    def get(self, species, genomic_datasets):
        """ Returns a list of subjects in the selected datasets  """
        args = filter_arguments.parse_args(request)

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in genomic_subject_filters.keys():
                raise BadRequest('Bad filter string %s' % args['filter'])

        if 'study_name' not in required_cols:
            required_cols = ['study_name'] + required_cols
        if 'dataset' not in required_cols:
            required_cols.append('dataset')
        if 'annotation_path' in required_cols:
            if 'annotation_method' not in required_cols:
                required_cols.append('annotation_method')
            if 'annotation_reference' not in required_cols:
                required_cols.append('annotation_reference')

        attribute_query = [genomic_subject_filters['id']['field']]        # the query requires the first field to be from Sample

        for col in required_cols:
            if col != 'id' and genomic_subject_filters[col]['field'] is not None:
                attribute_query.append(genomic_subject_filters[col]['field'])

        filter = json.loads(args['filter'])
        datasets = genomic_datasets.split(',')
        ret = find_genomic_subjects(attribute_query, species, datasets, filter)

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
                uniques['names_by_dataset'][s['dataset']].append(s['identifier'])

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
            sort_specs = [{'field': 'identifier', 'order': 'asc'}]

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
            'pages': ceil((total_size*1.0)/args['page_size'])
        }


def find_genomic_subjects(attribute_query, species, genomic_datasets, genomic_filters):
    results = []
    for dataset in genomic_datasets:
        db = get_genomic_db(species, dataset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        subject_query = db.session.query(*attribute_query)\
            .select_from(Subject)\
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
                        f['value'].append('')
                    filter_spec.append(f)
            except:
                raise BadRequest('Bad filter string %s')

        if len(filter_spec) > 0:
            subject_query = apply_filters(subject_query, filter_spec)

        if allele_filters is not None:
            subjects_with_alleles = db.session.query(Subject.identifier)\
                .join(SubjectSequence)\
                .join(Sequence).filter(Sequence.name.in_(allele_filters['value'])).all()
            subjects_with_alleles = [x[0] for x in subjects_with_alleles]
            subject_query = subject_query.filter(Subject.identifier.in_(subjects_with_alleles))

        subjects = subject_query.all()

        for s in subjects:
            r = s._asdict()
            for k in list(r.keys()):
                v = r[k]
                if isinstance(v, datetime):
                    r[k] = v.date().isoformat()
                elif k == 'annotation_path' and 'http' not in v:
                    r[k] = '/'.join([app.config['STATIC_LINK'], 'study_data/Genomic/samples', r[k]])
            r['dataset'] = dataset
            results.append(r)

    return results


def find_genomic_filter_params(species, genomic_datasets):
    genes = []

    for dset in genomic_datasets:
        db = get_genomic_db(species, dset)

        if db is None:
            raise BadRequest('Bad species or dataset name')

        g = db.session.query(Sequence.name)\
            .join(SubjectSequence, SubjectSequence.sequence_id == Sequence.id)\
            .join(Subject, SubjectSequence.subject_id == Subject.id)\
            .filter(Sequence.type.like('%REGION')) \
            .all()

        genes.extend(g)

    genes = sorted([gene[0] for gene in genes])
    gene_types = [gene[0:4] for gene in genes]
    gene_types = sorted(set(gene_types))

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


