# Services related to genomic sequences and features

from flask import request
from flask_restplus import Resource, reqparse, inputs
from api.restplus import api
from sqlalchemy import inspect

from app import db
from db.feature_db import Species, RefSeq, Feature, Sequence

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


filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('imgt', type=inputs.boolean)
filter_arguments.add_argument('novel', type=inputs.boolean)
filter_arguments.add_argument('full', type=inputs.boolean)
filter_arguments.add_argument('filter', type=str)
filter_arguments.add_argument('page_number', type=int)
filter_arguments.add_argument('page_size', type=int)

@ns.route('/sequences/<string:species>/<string:ref_seq>')
@api.response(404, 'Reference sequence not found.')
class SequencesAPI(Resource):
    @api.expect(filter_arguments, validate=True)
    def get(self, species, ref_seq):
        """ Returns nucleotide sequences from selected references. Use 'all' to wildcard species or ref_seq """

        args = filter_arguments.parse_args(request)
        ref_seqs = db.session.query(RefSeq)

        if species != 'all':
            ref_seqs = ref_seqs.join(Species).filter(Species.name == species)

        if ref_seq != 'all':
            ref_seqs = ref_seqs.filter(RefSeq.name == ref_seq)

        if not ref_seqs.count():
            return None, 404

        if ('page_size' in args and args['page_size'] <= 0) or ('page_number' in args and args['page_number'] < 0):
            return None, 404

        sequences = []

        for ref in ref_seqs.all():
            x_keys = {}
            if species == 'all':
                x_keys['species'] = ref.species.name
            if ref_seq == 'all':
                x_keys['ref_seq'] = ref.name
            for f in ref.features:
                for sequence in f.sequences:
                    if 'imgt' in args and not args['imgt']:
                        if not sequence.novel:
                            continue
                    if 'novel' in args and not args['novel']:
                        if sequence.novel:
                            continue
                    if 'full' in args and not args['full']:
                        if sequence.type not in ['V-REGION', 'D-REGION', 'J-REGION']:
                            continue
                    if 'filter' in args and args['filter']:
                        if args['filter'] not in sequence.name:
                            continue

                    seq = object_as_dict(sequence)
                    del seq['id']
                    del seq['species_id']
                    if len(x_keys):
                        seq.update(x_keys)
                    sequences.append(seq)

        if 'page_size' in args and 'page_number' in args:
            first = args['page_number'] * args['page_size']
            sequences = sequences[first : first + args['page_size']]

        return {'sequences': sequences}


@ns.route('/feature_pos/<string:species>/<string:ref_seq_name>/<string:feature_string>')
@api.response(404, 'Reference sequence not found.')
class FeaturePosAPI(Resource):
    def get(self, species, ref_seq_name, feature_string):
        """ Returns the position of the first feature matching the specified string """

        ref_seqs = db.session.query(RefSeq).join(Species).filter(Species.name == species.replace('_', ' ')).filter(RefSeq.name == ref_seq_name).all()       # fudge for space in species name

        if len(ref_seqs) != 1:
            return None, 404

        ref_seq = ref_seqs[0]
        q = db.session.query(Feature).join(RefSeq).filter(RefSeq.id == ref_seq.id)
        features = db.session.query(Feature).join(RefSeq).filter(RefSeq.id == ref_seq.id).filter(Feature.name.contains(feature_string)).all()

        start = 99999999
        end = 0

        for feature in features:
            start = min(start, feature.start)
            end = max(end, feature.end)

        return[{'chromosome': ref_seq_name, 'start': start, 'end': end}]

