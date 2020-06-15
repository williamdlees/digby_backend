# Services related to vdjbase repseq-based data sets

from flask import request
from flask_restplus import Resource, reqparse, fields, marshal, inputs
from api.restplus import api
from sqlalchemy import inspect, func
from math import ceil
import json
from sqlalchemy_filters import apply_filters
from werkzeug.exceptions import BadRequest
from datetime import datetime
import decimal


from app import vdjbase_dbs
from db.vdjbase_model import Sample, GenoDetection, Patient, SeqProtocol, Study, TissuePro, HaplotypesFile, SamplesHaplotype, Allele, AllelesSample

# Return SqlAlchemy row as a dict, using correct column names
def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


ns = api.namespace('repseq', description='Genes and annotations inferred from RepSeq data')

@ns.route('/species')
@api.response(404, 'No species available!')
class SpeciesApi(Resource):

    def get(self):

        """ Returns the list of species for which information is held """
        sp = list(vdjbase_dbs.keys())

        if sp:
            return sp
        else:
            return None, 404


@ns.route('/ref_seqs/<string:species>')
class DataSetAPI(Resource):

    def get(self, species):
        """ Returns the list of datasets available for the selected species """
        if species in vdjbase_dbs:
            return list(vdjbase_dbs[species].keys())
        else:
            return list()


valid_filters = {
    'name': {'model': 'Sample', 'field': Sample.name},
    'chain': {'model': 'Sample', 'field': Sample.chain},
    'row_reads': {'model': 'Sample', 'field': Sample.row_reads},
    'date': {'model': 'Sample', 'field': Sample.date},
    'samples_group': {'model': 'Sample', 'field': Sample.samples_group},

    'patient_name': {'model': 'Patient', 'field': Patient.name.label('patient_name'), 'fieldname': 'name'},
    'sex': {'model': 'Patient', 'field': Patient.sex},
    'ethnic': {'model': 'Patient', 'field': Patient.ethnic},
    'status': {'model': 'Patient', 'field': Patient.status},
    'cohort': {'model': 'Patient', 'field': Patient.cohort},
    'age': {'model': 'Patient', 'field': Patient.age},
    'name_in_paper': {'model': 'Patient', 'field': Patient.name_in_paper},
    'country': {'model': 'Patient', 'field': Patient.country},

    'study_name': {'model': 'Study', 'field': Study.name.label('study_name'), 'fieldname': 'name'},
    'institute': {'model': 'Study', 'field': Study.institute},
    'researcher': {'model': 'Study', 'field': Study.researcher},
    'num_subjects': {'model': 'Study', 'field': Study.num_subjects},
    'num_samples': {'model': 'Study', 'field': Study.num_samples},
    'study_reference': {'model': 'Study', 'field': Study.reference.label('study_reference'), 'fieldname': 'reference'},
    'study_contact': {'model': 'Study', 'field': Study.contact.label('study_contact'), 'fieldname': 'contact'},
    'study_acc_id': {'model': 'Study', 'field': Study.accession_id.label('study_acc_id'), 'fieldname': 'accession_id'},
    'study_acc_ref': {'model': 'Study', 'field': Study.accession_reference.label('study_acc_ref'), 'fieldname': 'accession_reference'},

    'tissue_name': {'model': 'TissuePro', 'field': TissuePro.name.label('tissue_name'), 'fieldname': 'name'},
    'species': {'model': 'TissuePro', 'field': TissuePro.species},
    'tissue': {'model': 'TissuePro', 'field': TissuePro.tissue},
    'combined_cell_type': {'model': 'TissuePro', 'field': TissuePro.combined_cell_type},
    'cell_type': {'model': 'TissuePro', 'field': TissuePro.cell_type},
    'sub_cell_type': {'model': 'TissuePro', 'field': TissuePro.sub_cell_type},
    'isotype': {'model': 'TissuePro', 'field': TissuePro.isotype},

    'sequencing_protocol': {'model': 'SeqProtocol', 'field': SeqProtocol.name.label('sequencing_protocol'), 'fieldname': 'name'},
    'umi': {'model': 'SeqProtocol', 'field': SeqProtocol.umi},
    'sequencing_length': {'model': 'SeqProtocol', 'field': SeqProtocol.sequencing_length},
    'primer_3_loc': {'model': 'SeqProtocol', 'field': SeqProtocol.primers_3_location.label('primer_3_loc'), 'fieldname': 'primers_3_location'},
    'primer_5_loc': {'model': 'SeqProtocol', 'field': SeqProtocol.primers_5_location.label('primer_5_loc'), 'fieldname': 'primers_5_location'},
    'seq_platform': {'model': 'SeqProtocol', 'field': SeqProtocol.sequencing_platform.label('seq_platform'), 'fieldname': 'sequencing_platform'},
    'helix': {'model': 'SeqProtocol', 'field': SeqProtocol.helix},

    'pipeline_name': {'model': 'GenoDetection', 'field': GenoDetection.name.label('pipeline_name'), 'fieldname': 'name'},
    'prepro_tool': {'model': 'GenoDetection', 'field': GenoDetection.prepro_tool},
    'aligner_tool': {'model': 'GenoDetection', 'field': GenoDetection.aligner_tool},
    'aligner_ver': {'model': 'GenoDetection', 'field': GenoDetection.aligner_ver},
    'aligner_reference': {'model': 'GenoDetection', 'field': GenoDetection.aligner_reference},
    'geno_tool': {'model': 'GenoDetection', 'field': GenoDetection.geno_tool},
    'geno_ver': {'model': 'GenoDetection', 'field': GenoDetection.geno_ver},
    'haplotype_tool': {'model': 'GenoDetection', 'field': GenoDetection.haplotype_tool},
    'haplotype_ver': {'model': 'GenoDetection', 'field': GenoDetection.haplotype_ver},
    'single_assignment': {'model': 'GenoDetection', 'field': GenoDetection.single_assignment},
    'detection': {'model': 'GenoDetection', 'field': GenoDetection.detection},

    'allele': {'model': None, 'field': None},

    'haplotypes': {'model': None, 'field': None},
    'genotypes': {'model': None, 'field': None}
}



@ns.route('/sample_info/<string:species>/<string:dataset>/<string:sample>')
class SampleInfoApi(Resource):
    def get(self, species, dataset, sample):
        """ Returns information on the selected sample """

        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        session = vdjbase_dbs[species][dataset].session
        attribute_query = []

        for col in valid_filters.keys():
            if valid_filters[col]['field'] is not None:
                attribute_query.append(valid_filters[col]['field'])

        info = session.query(*attribute_query).join(GenoDetection).join(Patient).join(SeqProtocol).join(TissuePro).join(Study, Sample.study_id == Study.id).filter(Sample.name==sample).one_or_none()

        if info:
            info = info._asdict()
            info['date'] = info['date'].isoformat()
            haplotypes = session.query(HaplotypesFile.by_gene_s).join(SamplesHaplotype).join(Sample).filter(Sample.name==sample).order_by(HaplotypesFile.by_gene_s).all()
            info['haplotypes'] = [(h[0]) for h in haplotypes]

        return info



filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('page_number', type=int)
filter_arguments.add_argument('page_size', type=int)
filter_arguments.add_argument('filter', type=str)
filter_arguments.add_argument('sort_by', type=str)
filter_arguments.add_argument('cols', type=str)


@ns.route('/samples/<string:species>/<string:dataset>')
class SamplesApi(Resource):
    @api.expect(filter_arguments, validate=True)
    def get(self, species, dataset):
        """ Returns the list of samples in the selected dataset """

        if species not in vdjbase_dbs or set(dataset.split(',')).difference(set(vdjbase_dbs[species])):
            return list()

        args = filter_arguments.parse_args(request)
        if (args['page_size'] and args['page_size'] <= 0) or (args['page_number'] and args['page_number'] < 0):
            return list()

        if args['cols'] is None or not args['cols']:
            return list()

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in valid_filters.keys():
                print('bad column in request: %s' % col)
                return list(), 404

        hap_filters = None
        allele_filters = None

        filter_spec = []
        if args['filter']:
            for f in json.loads(args['filter']):
                try:
                    if f['field'] == 'haplotypes':
                        hap_filters = f
                    elif f['field'] == 'allele':
                        allele_filters = f
                    else:
                        f['model'] = valid_filters[f['field']]['model']
                        if 'fieldname' in valid_filters[f['field']]:
                            f['field'] = valid_filters[f['field']]['fieldname']
                        if '(blank)' in f['value']:
                            f['value'].append('')
                        filter_spec.append(f)
                except:
                    raise BadRequest('Bad filter string %s' % args['filter'])

        ret = []

        for dset in dataset.split(','):
            session = vdjbase_dbs[species][dset].session

            attribute_query = []

            for col in required_cols:
                if valid_filters[col]['field'] is not None:
                    attribute_query.append(valid_filters[col]['field'])

            query = session.query(*attribute_query).join(GenoDetection).join(Patient).join(SeqProtocol).join(TissuePro).join(Study, Sample.study_id == Study.id)
            query = apply_filters(query, filter_spec)

            if hap_filters:
                hap_samples = session.query(Sample.name.distinct()).join(SamplesHaplotype).join(HaplotypesFile).filter(HaplotypesFile.by_gene_s.in_(hap_filters['value']))
                query = query.filter(Sample.name.in_(hap_samples))

            if allele_filters:
                allele_samples = session.query(Sample.name.distinct()).join(AllelesSample, Sample.id == AllelesSample.sample_id).join(Allele, Allele.id == AllelesSample.allele_id).filter(Allele.name.in_(allele_filters['value'])).all()
                if allele_samples is None:
                    allele_samples = []
                query = query.filter(Sample.name.in_([s[0] for s in allele_samples]))


            res = query.all()

            for r in res:
                s = r._asdict()
                for k, v in s.items():
                    if isinstance(v, datetime):
                        s[k] = v.date().isoformat()
                s['dataset'] = dset
                ret.append(s)

        total_size = len(ret)

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in ret:
            for f in required_cols:
                if valid_filters[f]['field'] is not None and 'no_uniques' not in valid_filters[f]:
                    el = s[f]
                    if isinstance(el, datetime):
                        s[f] = el.date().isoformat()
                    elif isinstance(el, str) and len(el) == 0:
                        s[f] = '(blank)'
                    if s[f] not in uniques[f]:
                        uniques[f].append(s[f])

        for f in required_cols:
            try:
                uniques[f].sort(key=lambda x: (x is None or x == '', x))
            except:
                pass

        if 'haplotypes' in required_cols:
            haplotypes = session.query(HaplotypesFile.by_gene_s).distinct().order_by(HaplotypesFile.by_gene_s).all()
            uniques['haplotypes'] = [(h[0]) for h in haplotypes]

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        for spec in sort_specs:
            if spec['field'] in valid_filters.keys():
                ret = sorted(ret, key=lambda x : ((x[spec['field']] is None or x[spec['field']] == ''),  x[spec['field']]), reverse=(spec['order'] == 'desc'))

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first : first + args['page_size']]

        if 'haplotypes' in required_cols:                                   # BUG - this won't work with multiple datasets
            haplotypes = session.query(Sample.name, func.group_concat(HaplotypesFile.by_gene_s))

            for r in ret:
                h = haplotypes.filter(Sample.name == r['name']).join(SamplesHaplotype).join(HaplotypesFile).one_or_none()
                if h:
                    r['haplotypes'] = h[1]
                else:
                    r['haplotypes'] = ''

        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size'])
        }

valid_sequence_cols = {
    'name': {'model': 'Allele', 'field': Allele.name},
    'seq': {'model': 'Allele', 'field': Allele.seq, 'no_uniques': True},
    'seq_len': {'model': 'Allele', 'field': Allele.seq_len},
    'similar': {'model': 'Allele', 'field': Allele.similar},
    'appears': {'model': 'Allele', 'field': Allele.appears},
    'is_single_allele': {'model': 'Allele', 'field': Allele.is_single_allele},
    'low_confidence': {'model': 'Allele', 'field': Allele.low_confidence},
    'novel': {'model': 'Allele', 'field': Allele.novel},
    'max_kdiff': {'model': 'Allele', 'field': Allele.max_kdiff},
}

@ns.route('/sequences/<string:species>/<string:dataset>')
class SequencesApi(Resource):
    @api.expect(filter_arguments, validate=True)
    def get(self, species, dataset):
        """ Returns the list of sequences in the selected dataset """

        if species not in vdjbase_dbs or set(dataset.split(',')).difference(set(vdjbase_dbs[species])):
            return list()

        args = filter_arguments.parse_args(request)
        if (args['page_size'] and args['page_size'] <= 0) or (args['page_number'] and args['page_number'] < 0):
            return list()

        if args['cols'] is None or not args['cols']:
            return list()

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in valid_sequence_cols.keys():
                print('bad column in request: %s' % col)
                return list(), 404

        filter_spec = []
        if args['filter']:
            for f in json.loads(args['filter']):
                try:
                    if f['field'] != 'haplotypes':
                        f['model'] = valid_sequence_cols[f['field']]['model']
                        if 'fieldname' in valid_sequence_cols[f['field']]:
                            f['field'] = valid_sequence_cols[f['field']]['fieldname']
                        if '(blank)' in f['value']:
                            f['value'].append('')
                        filter_spec.append(f)
                    else:
                        hap_filters = f
                except:
                    raise BadRequest('Bad filter string %s' % args['filter'])

        ret = []

        for dset in dataset.split(','):
            session = vdjbase_dbs[species][dset].session

            attribute_query = []

            for col in required_cols:
                if valid_sequence_cols[col]['field'] is not None:
                    attribute_query.append(valid_sequence_cols[col]['field'])

            query = session.query(*attribute_query)
            query = apply_filters(query, filter_spec)
            res = query.all()

            for r in res:
                s = r._asdict()
                for k, v in s.items():
                    if isinstance(v, datetime):
                        s[k] = v.date().isoformat()
                    elif isinstance(v, decimal.Decimal):
                        s[k] = '%0.2f' % v
                s['dataset'] = dset
                ret.append(s)

        total_size = len(ret)

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in ret:
            for f in required_cols:
                if valid_sequence_cols[f]['field'] is not None and 'no_uniques' not in valid_sequence_cols[f]:
                    el = s[f]
                    if isinstance(el, datetime):
                        s[f] = el.date().isoformat()
                    elif isinstance(el, decimal.Decimal):
                        s[f] = '%0.2f' % el
                    elif isinstance(el, str) and len(el) == 0:
                        s[f] = '(blank)'
                    if s[f] not in uniques[f]:
                        uniques[f].append(s[f])

        for f in required_cols:
            try:
                uniques[f].sort(key=lambda x: (x is None or x == '', x))
            except:
                pass

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        for spec in sort_specs:
            if spec['field'] in valid_sequence_cols.keys():
                ret = sorted(ret, key=lambda x : ((x[spec['field']] is None or x[spec['field']] == ''),  x[spec['field']]), reverse=(spec['order'] == 'desc'))

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
