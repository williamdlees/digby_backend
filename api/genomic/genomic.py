# Services related to genomic sequences and features

from flask import request
from flask_restplus import Resource, reqparse, inputs
from api.restplus import api
from sqlalchemy import inspect, or_, func, distinct
from math import ceil
from werkzeug.exceptions import BadRequest
from app import db
from db.feature_db import Species, RefSeq, Feature, Sequence, SequenceFeature, Sample, Study, SampleSequence
import json
from datetime import datetime
from sqlalchemy_filters import apply_filters
from decimal import Decimal
from app import app

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

        if sp:
            return [row.name for row in sp]
        else:
            return None, 404


@ns.route('/ref_seqs/<string:species>')
@api.response(404, 'No reference sequences available for that species.')
class RefSeqAPI(Resource):

    def get(self, species):
        """ Returns the list of annotated reference sequences for the selected species """
        sp = db.session.query(Species).filter_by(name=species).one_or_none()

        if sp:
            return [{'ref_seq': row.name, 'locus': row.locus} for row in sp.refseqs]
        else:
            return None, 404

@ns.route('/sample_info/<string:species>/<string:study_name>/<string:sample>')
class SampleInfoApi(Resource):
    def get(self, species, study_name, sample):
        """ Returns information on the selected sample """

        sample = db.session.query(Sample)\
            .join(Species)\
            .join(Study, Study.id == Sample.study_id)\
            .join(RefSeq, RefSeq.id == Sample.ref_seq_id)\
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
            .join(Species)\
            .join(Study, Study.id == Sample.study_id)\
            .join(RefSeq, RefSeq.id == Sample.ref_seq_id)\
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
    'name': {'model': 'Sequence', 'field': Sequence.name},
    'imgt_name': {'model': 'Sequence', 'field': Sequence.imgt_name},
    'type': {'model': 'Sequence', 'field': Sequence.type},
    'novel': {'model': 'Sequence', 'field': Sequence.novel},
    'deleted': {'model': 'Sequence', 'field': Sequence.deleted},
    'sequence': {'model': 'Sequence', 'field': Sequence.sequence, 'no_uniques': True},
    'gapped_sequence': {'model': 'Sequence', 'field': Sequence.gapped_sequence, 'no_uniques': True},

    'feature': {'model': 'Feature', 'field': Feature.feature},
    'start': {'model': 'Feature', 'field': Feature.start},
    'end': {'model': 'Feature', 'field': Feature.end},
    'score': {'model': 'Feature', 'field': Feature.score},
    'strand': {'model': 'Feature', 'field': Feature.strand},
    'frame': {'model': 'Feature', 'field': Feature.frame},

    'refseq_name': {'model': 'RefSeq', 'field': RefSeq.name.label('refseq_name'), 'fieldname': 'name'},

    'sample_count': {'field': func.count(Sample.name).label('sample_count'), 'fieldname': 'sample_count'},
    'appearances': {'field': func.sum(SampleSequence.chromo_count).label('appearances'), 'fieldname': 'appearances'},

    'sample_id': {'model': None, 'fieldname': 'sample_id'},
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

        refs = db.session.query(RefSeq.id).join(Species).filter(Species.name == species).filter(RefSeq.name.in_(genomic_datasets.split(','))).all()

        if not refs:
            raise BadRequest('Bad species name or reference set name %s' % species)

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in genomic_sequence_filters.keys():
                raise BadRequest('Bad column string %s' % args['cols'])

        attribute_query = [genomic_sequence_filters['name']['field']]        # the query requires the first field to be from Sequence

        for col in required_cols:
            if col != 'name' and genomic_sequence_filters[col]['field'] is not None:
                attribute_query.append(genomic_sequence_filters[col]['field'])

        seq_query = db.session.query(*attribute_query)\
            .join(SequenceFeature, SequenceFeature.sequence_id == Sequence.id)\
            .join(Feature)\
            .join(RefSeq)\
            .join(SampleSequence, SampleSequence.sequence_id == Sequence.id)\
            .join(Sample)\
            .filter(RefSeq.id.in_(refs))\
            .group_by(Sequence.name)

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
                        if '(blank)' in f['value']:
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
            alleles_with_samples = db.session.query(Sequence.name)\
                .join(SampleSequence)\
                .join(Sample)\
                .filter(Sample.id.in_(sample_id_filter['value'])).all()
            seq_query = seq_query.filter(Sequence.name.in_(alleles_with_samples))

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
                    if el not in uniques[f]:
                        uniques[f].append(el)

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

    'study_name': {'model': 'Study', 'field': Study.name.label('study_name'), 'fieldname': 'study_name'},
    'institute': {'model': 'Study', 'field': Study.institute},
    'researcher': {'model': 'Study', 'field': Study.researcher},
    'reference': {'model': 'Study', 'field': Study.reference},
    'contact': {'model': 'Study', 'field': Study.contact},
    'accession_id': {'model': 'Study', 'field': Study.accession_id},
    'accession_reference': {'model': 'Study', 'field': Study.accession_reference},

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

        attribute_query = [genomic_sample_filters['id']['field']]        # the query requires the first field to be from Sample

        for col in required_cols:
            if col != 'id' and genomic_sample_filters[col]['field'] is not None:
                attribute_query.append(genomic_sample_filters[col]['field'])

        samples = find_genomic_samples(attribute_query, species, genomic_datasets.split(','), json.loads(args['filter']))

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in samples:
            for f in required_cols:
                if genomic_sample_filters[f]['field'] is not None and 'no_uniques' not in genomic_sample_filters[f]:
                    el = getattr(s, f)
                    if isinstance(el, datetime):
                        el = el.date().isoformat()
                    elif isinstance(el, str) and len(el) == 0:
                        el = '(blank)'
                    elif f == 'report' and 'http' not in el:
                        el = app.config['STATIC_LINK'] + el.replace('\\', '/').replace(' ', '%20')
                    if el not in uniques[f]:
                        uniques[f].append(el)

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
        ref = db.session.query(RefSeq.id)\
            .join(Species)\
            .filter(Species.id == sp.id)\
            .filter(RefSeq.name == dset)
        ref = ref.one_or_none()
        if ref is None:
            return 'No such genomic dataset %s' % dset, '404'
        ref_ids.append(ref.id)

    sample_query = db.session.query(*attribute_query)\
        .join(Study)\
        .filter(Sample.ref_seq_id.in_(ref_ids))

    allele_filters = None

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

