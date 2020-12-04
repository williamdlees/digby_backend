# Services related to genomic sequences and features

from flask import request
from flask_restx import Resource, reqparse, inputs
from api.restx import api
from sqlalchemy import inspect, or_, func, distinct
from math import ceil
from werkzeug.exceptions import BadRequest
from app import sql_db as db
from db.feature_db import Species, RefSeq, Feature, Sequence, SequenceFeature, Sample, Study, SampleSequence, DataSet
import json
from datetime import datetime
from sqlalchemy_filters import apply_filters
from decimal import Decimal
from app import app
from api.vdjbase.vdjbase import get_vdjbase_species

# Return SqlAlchemy row as a dict, using correct column names
def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


ns = api.namespace('genomic', description='Genomic data and annotations')

@ns.route('/species')
@api.response(404, 'No species available!')
class SpeciesApi(Resource):

    def get(self):
        """ Returns the list of species for which information is held """
        sp = db.session.query(Species).all()
        genomic_sp = [row.name for row in sp] if sp else []
        vdjbase_sp = get_vdjbase_species()
        return list(set(genomic_sp) | set(vdjbase_sp))


@ns.route('/assemblies/<string:species>/<string:data_sets>')
class AssemblyAPI(Resource):

    def get(self, species, data_sets):
        """ Returns the list of annotated assemblies for the selected species and datasets """
        data_sets = data_sets.split(',')

        assemblies = db.session.query(distinct(RefSeq.name))\
            .join(Sample)\
            .join(DataSet)\
            .join(Species)\
            .filter(Species.name == species)\
            .filter(DataSet.name.in_(data_sets)).all()

        return [{'assembly': row[0]} for row in assemblies]


@ns.route('/data_sets/<string:species>')
@api.response(404, 'No data sets available for that species.')
class DataSetAPI(Resource):

    def get(self, species):
        """ Returns the list of data sets for the selected species """
        sp = db.session.query(Species).filter_by(name=species).one_or_none()

        if sp:
            return [{'dataset': row.name, 'locus': row.locus} for row in sp.datasets]
        else:
            return []

@ns.route('/sample_info/<string:species>/<string:study_name>/<string:sample>')
class SampleInfoApi(Resource):
    def get(self, species, study_name, sample):
        """ Returns information on the selected sample """

        sample = db.session.query(Sample)\
            .join(Species)\
            .join(Study, Study.id == Sample.study_id)\
            .filter(Species.name == species)\
            .filter(Study.name == study_name)\
            .filter(Sample.name == sample).one_or_none()

        if sample is None:
            raise BadRequest('Bad species name, study ID or sample name')

        attribute_query = []

        for col in genomic_sample_filters.keys():
            if genomic_sample_filters[col]['field'] is not None:
                attribute_query.append(genomic_sample_filters[col]['field'])

        info = db.session.query(*attribute_query)\
            .join(Species, Species.id == Sample.species_id)\
            .join(Study, Study.id == Sample.study_id)\
            .join(RefSeq, RefSeq.id == Sample.ref_seq_id)\
            .join(DataSet, DataSet.id == Sample.data_set_id)\
            .filter(Species.name == species)\
            .filter(Study.name == study_name)\
            .filter(Sample.id == sample.id).one_or_none()

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
    'name': {'model': 'Sequence', 'field': Sequence.name, 'sort': 'name'},
    'imgt_name': {'model': 'Sequence', 'field': Sequence.imgt_name},
    'type': {'model': 'Sequence', 'field': Sequence.type},
    'novel': {'model': 'Sequence', 'field': Sequence.novel},
    'deleted': {'model': 'Sequence', 'field': Sequence.deleted},
    'functional': {'model': 'Sequence', 'field': Sequence.functional},
    'sequence': {'model': 'Sequence', 'field': Sequence.sequence, 'no_uniques': True},
    'gapped_sequence': {'model': 'Sequence', 'field': Sequence.gapped_sequence, 'no_uniques': True},

    'appearances': {'field': func.count(Sample.id.distinct()).label('appearances'), 'fieldname': 'appearances', 'sort': 'numeric'},
    'sample_id': {'model': None, 'fieldname': 'sample_id'},
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
    @api.expect(filter_arguments, validate=True)
    def get(self, species, genomic_datasets):
        """ Returns nucleotide sequences from selected reference or multiple references (separate multiple reference names with ',')  """
        args = filter_arguments.parse_args(request)

        sp = db.session.query(Species).filter(Species.name == species).one_or_none()
        if sp is None:
            raise BadRequest('Bad species name %s' % species)

        ds = db.session.query(DataSet.id)\
            .filter(DataSet.species_id == sp.id)\
            .filter(DataSet.name.in_(genomic_datasets.split(',')))\
            .all()
        if ds is None:
            raise BadRequest('Bad data set name in %s' % genomic_datasets)

        sequence_id_query = db.session.query(Sequence.id.distinct())\
            .join(SampleSequence, Sequence.id == SampleSequence.sequence_id)\
            .join(Sample, Sample.id == SampleSequence.sample_id)\
            .join(DataSet)\
            .join(Species)\
            .filter(Species.name == species)\
            .filter(DataSet.id.in_(ds))

        sequence_ids = sequence_id_query.all()

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in genomic_sequence_filters.keys():
                raise BadRequest('Bad column string %s' % args['cols'])

        if 'sequence' in required_cols and 'gapped_sequence' not in required_cols:
            required_cols.append('gapped_sequence')

        attribute_query = [genomic_sequence_filters['name']['field']]        # the query requires the first field to be from Sequence

        for col in required_cols:
            if col != 'name' and genomic_sequence_filters[col]['field'] is not None:
                attribute_query.append(genomic_sequence_filters[col]['field'])

        seq_query = db.session.query(*attribute_query)\
            .filter(Sequence.id.in_(sequence_ids))\
            .join(SampleSequence, Sequence.id == SampleSequence.sequence_id)\
            .join(Sample, Sample.id == SampleSequence.sample_id)\
            .group_by(Sequence.id)


        filter_spec = []
        sample_count_filters = []
        sample_id_filter = None
        appearances_filters = []
        if args['filter']:
            for f in json.loads(args['filter']):
                try:
                    if 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'sample_count':
                        sample_count_filters.append(f)
                    elif 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'sample_id':
                        sample_id_filter = f
                    elif 'fieldname' in genomic_sequence_filters[f['field']] and genomic_sequence_filters[f['field']]['fieldname'] == 'appearances':
                        appearances_filters.append(f)
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
                    raise BadRequest('Bad filter string %s' % args['filter'])

        seq_query = apply_filters(seq_query, filter_spec)

        for f in sample_count_filters:
            if f['op'] in OPERATORS:
                seq_query = seq_query.having(OPERATORS[f['op']](func.count(Sample.name), f['value']))

        for f in appearances_filters:
            if f['op'] in OPERATORS:
                seq_query = seq_query.having(OPERATORS[f['op']](func.sum(SampleSequence.chromo_count), f['value']))

        if sample_id_filter is not None:
            filtered_sample_ids = []

            for dataset, names in sample_id_filter['value'].items():
                sample_ids = db.session.query(Sample.id.distinct())\
                    .join(DataSet)\
                    .filter(DataSet.name == dataset)\
                    .filter(Sample.name.in_(names)).all()
                filtered_sample_ids.extend(sample_ids)

            filtered_sequence_ids = sequence_id_query.filter(Sample.id.in_(filtered_sample_ids)).all()
            seq_query = seq_query.filter(Sequence.id.in_(filtered_sequence_ids))

        seqs = seq_query.all()

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in seqs:
            for f in required_cols:
                if genomic_sequence_filters[f]['field'] is not None and 'no_uniques' not in genomic_sequence_filters[f]:
                    el = getattr(s, f)
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

        # something of a sort order for names

        def allele_sort_key(x):
            if x is None or x == '' or 'IGH' not in x:
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
                elif 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'name':
                    uniques[f].sort(key=allele_sort_key)
                else:
                    uniques[f].sort(key=lambda x: (x is None or x == '', x))
            except:
                pass

        ret = []
        for r in seqs:
            s = r._asdict()
            for k, v in s.items():
                if isinstance(v, datetime):
                    s[k] = v.date().isoformat()
                elif isinstance(v, Decimal):
                    s[k] = int(v)

            ret.append(s)

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        for spec in sort_specs:
            if spec['field'] in genomic_sequence_filters.keys():
                f = spec['field']
                if 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'numeric':
                    ret = sorted(ret, key=lambda x: num_sort_key(x[spec['field']]), reverse=(spec['order'] == 'desc'))
                elif 'sort' in genomic_sequence_filters[f] and genomic_sequence_filters[f]['sort'] == 'name':
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


@ns.route('/feature_pos/<string:species>/<string:ref_seq_name>/<string:feature_string>')
@api.response(404, 'Reference sequence not found.')
class FeaturePosAPI(Resource):
    def get(self, species, ref_seq_name, feature_string):
        """ Returns the position of the first feature matching the specified string """

        ref_seqs = db.session.query(RefSeq).join(Species).filter(Species.name == species.replace('_', ' ')).filter(RefSeq.name == ref_seq_name).all()       # fudge for space in species name

        if len(ref_seqs) != 1:
            return None, 404

        ref_seq = ref_seqs[0]
        features = db.session.query(Feature).join(RefSeq).join(Species).filter(Species.name == species.replace('_', ' ')).filter(RefSeq.id == ref_seq.id).filter(Feature.name.contains(feature_string)).all()

        if features:
            start = 99999999
            end = 0

            for feature in features:
                start = min(start, feature.start)
                end = max(end, feature.end)

            return[{'chromosome': ref_seq_name, 'start': start, 'end': end}]
        else:
            return []


genomic_sample_filters = {
    'name': {'model': 'Sample', 'field': Sample.name},
    'id': {'model': 'Sample', 'field': Sample.id},
    'type': {'model': 'Sample', 'field': Sample.type},
    'date': {'model': 'Sample', 'field': Sample.date},
    'report': {'model': 'Sample', 'field': Sample.report_link.label('report'), 'fieldname': 'report'},
    'sample_description': {'model': 'Sample', 'field': Sample.description.label('sample_description'), 'fieldname': 'sample_description'},

    'study_name': {'model': 'Study', 'field': Study.name.label('study_name'), 'fieldname': 'study_name'},
    'institute': {'model': 'Study', 'field': Study.institute},
    'researcher': {'model': 'Study', 'field': Study.researcher},
    'reference': {'model': 'Study', 'field': Study.reference},
    'contact': {'model': 'Study', 'field': Study.contact},
    'study_description': {'model': 'Study', 'field': Study.description.label('study_description'), 'fieldname': 'study_description'},

    'assembly_id': {'model': 'RefSeq', 'field': RefSeq.name.label('assembly_id'), 'fieldname': 'assembly_id'},
    'assembly_reference': {'model': 'RefSeq', 'field': RefSeq.reference.label('assembly_reference'), 'fieldname': 'assembly_reference'},
    'chromosome': {'model': 'RefSeq', 'field': RefSeq.chromosome},
    'assembly_start': {'model': 'RefSeq', 'field': RefSeq.start.label('assembly_start'), 'fieldname': 'assembly_start', 'sort': 'numeric'},
    'assembly_end': {'model': 'RefSeq', 'field': RefSeq.end.label('assembly_end'), 'fieldname': 'assembly_end', 'sort': 'numeric'},

    'dataset': {'model': 'Dataset', 'field': DataSet.name.label('dataset'), 'fieldname': 'dataset'},

    'allele': {'model': None, 'field': None},

}


@ns.route('/samples/<string:species>/<string:genomic_datasets>')
@api.response(404, 'Reference sequence not found.')
class SamplesAPI(Resource):
    @api.expect(filter_arguments, validate=True)
    def get(self, species, genomic_datasets):
        """ Returns a list of samples that provide results against the selected reference or multiple references (separate multiple reference names with ',')  """
        args = filter_arguments.parse_args(request)

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in genomic_sample_filters.keys():
                raise BadRequest('Bad filter string %s' % args['filter'])

        if 'study_name' not in required_cols:
            required_cols = ['study_name'] + required_cols
        if 'dataset' not in required_cols:
            required_cols.append('dataset')
        if 'type' not in required_cols:
            required_cols.append('type')

        attribute_query = [genomic_sample_filters['id']['field']]        # the query requires the first field to be from Sample

        for col in required_cols:
            if col != 'id' and genomic_sample_filters[col]['field'] is not None:
                attribute_query.append(genomic_sample_filters[col]['field'])

        filter = json.loads(args['filter'])
        datasets = genomic_datasets.split(',')
        samples = find_genomic_samples(attribute_query, species, datasets, filter)

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

        for f in required_cols:
            try:
                if 'sort' in genomic_sample_filters[f] and genomic_sample_filters[f]['sort'] == 'numeric':
                    uniques[f].sort(key=num_sort_key)
                else:
                    uniques[f].sort(key=lambda x: (x is None or x == '', x))
            except:
                pass

        for s in samples:
            for f in required_cols:
                if genomic_sample_filters[f]['field'] is not None and 'no_uniques' not in genomic_sample_filters[f]:
                    el = getattr(s, f)
                    if isinstance(el, datetime):
                        el = el.date().isoformat()
                    elif isinstance(el, str) and len(el) == 0:
                        el = '(blank)'
                    elif f == 'report':
                        sample_type = getattr(s, 'type')
                        if type == 'Igenotyper assembly':
                            if 'http' not in el:
                                el = app.config['STATIC_LINK'] + el.replace('\\', '/').replace(' ', '%20')
                        elif type == 'VDJbase_assembly':
                           el = app.config['STATIC_LINK'] + el.replace('\\', '/').replace(' ', '%20')
                    if el not in uniques[f]:
                        uniques[f].append(el)
            if filter_applied:
                uniques['names_by_dataset'][s.dataset].append(s.name)

        ret = []
        for r in samples:
            s = r._asdict()
            for k in list(s.keys()):
                v = s[k]
                if isinstance(v, datetime):
                    s[k] = v.date().isoformat()
                elif k == 'report' and 'http' not in v:
                    s[k] = app.config['STATIC_LINK'] + s[k]
                    s['hap_report'] = s[k].replace('IGenotyper_report.html','genes_assigned_to_alleles.txt')
            ret.append(s)

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        for spec in sort_specs:
            if spec['field'] in genomic_sample_filters.keys():
                ret = sorted(ret, key=lambda x : (x[spec['field']] is None,  x[spec['field']]), reverse=(spec['order'] == 'desc'))

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

def find_genomic_samples(attribute_query, species, genomic_datasets, genomic_filters):
    sp = db.session.query(Species.id).filter(Species.name == species).one_or_none()

    if sp is None:
        raise BadRequest('No such species')

    ref_ids = []

    for dset in genomic_datasets:
        ref = db.session.query(DataSet.id)\
            .join(Species)\
            .filter(Species.id == sp.id)\
            .filter(DataSet.name == dset)
        ref = ref.one_or_none()
        if ref is None:
            raise BadRequest('No such genomic dataset %s' % dset)
        ref_ids.append(ref.id)

    sample_query = db.session.query(*attribute_query)\
        .join(Study)\
        .join(RefSeq, Sample.ref_seq_id == RefSeq.id)\
        .join(DataSet, Sample.data_set_id == DataSet.id)\
        .filter(Sample.data_set_id.in_(ref_ids))

    allele_filters = None
    dataset_filters = []

    filter_spec = []
    for f in genomic_filters:
        try:
            if f['field'] == 'allele':
                allele_filters = f
            else:
                f['model'] = genomic_sample_filters[f['field']]['model']
                if 'fieldname' in genomic_sample_filters[f['field']]:
                    f['field'] = genomic_sample_filters[f['field']]['fieldname']
                if '(blank)' in f['value']:
                    f['value'].append('')
                filter_spec.append(f)
        except:
            raise BadRequest('Bad filter string %s')

    sample_query = apply_filters(sample_query, filter_spec)

    if allele_filters is not None:
        samples_with_alleles = db.session.query(Sample.name)\
            .join(SampleSequence)\
            .join(Sequence).filter(Sequence.name.in_(allele_filters['value'])).all()
        sample_query = sample_query.filter(Sample.name.in_(samples_with_alleles))

    return sample_query.all()

def find_genomic_filter_params(species, genomic_datasets):
    sp = db.session.query(Species.id).filter(Species.name == species).one_or_none()

    if sp is None:
        return []

    for dset in genomic_datasets:
        ref = db.session.query(DataSet.id)\
            .join(Species)\
            .filter(Species.id == sp.id)\
            .filter(DataSet.name == dset)
        ref = ref.one_or_none()
        if ref is None:
            return 'No such genomic dataset %s' % dset, '404'

    genes = db.session.query(Sequence.name)\
        .join(SampleSequence, SampleSequence.sequence_id == Sequence.id)\
        .join(Sample, SampleSequence.sample_id == Sample.id)\
        .join(DataSet)\
        .filter(Sequence.type.like('%REGION')) \
        .filter(DataSet.name.in_(genomic_datasets)).all()

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
