# Services related to vdjbase repseq-based data sets

from flask import request
from flask_restplus import Resource, reqparse, fields, marshal, inputs
from api.restplus import api
from sqlalchemy import inspect
from math import ceil

from app import vdjbase_dbs
from db.vdjbase_model import Sample, GenoDetection, Patient, SeqProtocol, Study, TissuePro

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
@api.response(404, 'No reference sequences available for that species.')
class DataSetAPI(Resource):

    def get(self, species):
        """ Returns the list of datasets available for the selected species """
        if species in vdjbase_dbs:
            return list(vdjbase_dbs[species].keys())
        else:
            return None, 404

sample_model = api.model('Model', {
    'name': fields.String,
    'row_reads': fields.Integer,
    'genotype': fields.String,
    'genotype_graph': fields.String,
    'date': fields.DateTime(dt_format='rfc822'),
    'samples_group': fields.Integer,
})

filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('page_number', type=int)
filter_arguments.add_argument('page_size', type=int)


@ns.route('/samples/<string:species>/<string:dataset>')
@api.response(404, 'No such dataset.')
class SamplesApi(Resource):
    @api.expect(filter_arguments, validate=True)
    def get(self, species, dataset):
        """ Returns the list of samples in the selected dataset """

        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404


        samples = vdjbase_dbs[species][dataset].session.query(Sample).join(GenoDetection).join(Patient).join(SeqProtocol).join(TissuePro).join(Study, Sample.study_id == Study.id)

        args = filter_arguments.parse_args(request)
        if (args['page_size'] and args['page_size'] <= 0) or (args['page_number'] and args['page_number'] < 0):
            return None, 404

        if args['page_size'] and args['page_number']:
            first = args['page_number'] * args['page_size']
            samples = samples.slice(first, first + args['page_size'])
        else:
            samples = samples.all()

        ret = []

        for sample in samples:
            sample_dict = marshal(object_as_dict(sample), sample_model)
            sample_dict['geno_detection'] = object_as_dict(sample.geno_detection)
            sample_dict['patient'] = object_as_dict(sample.patient)
            sample_dict['seq_protocol'] = object_as_dict(sample.seq_protocol)
            sample_dict['study'] = object_as_dict(sample.study)
            sample_dict['tissue_pro'] = object_as_dict(sample.tissue_pro)
            ret.append(sample_dict)

        return ret
