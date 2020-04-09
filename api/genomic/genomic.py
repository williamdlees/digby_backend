# Services related to genomic sequences and features

from flask import request
from flask_restplus import Resource, reqparse
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


@ns.route('/jbrowse/<string:species>/features/<string:ref_seq>')
class JbrowseFeatureAPI(Resource):
    def add_subfeatures(self, feature, fq):
        if 'uniqueID' in feature:
            subs = fq.filter(Feature.parent_id == int(feature['uniqueID'])).all()

            if subs:
                feature['subfeatures'] = []
                for sub in subs:
                    en_sub = enumerate_feature(sub)
                    self.add_subfeatures(en_sub, fq)
                    feature['subfeatures'].append(en_sub)

    @api.expect(range_arguments, validate=True)
    @api.response(404, 'No features.')
    def get(self, species, ref_seq):
        """ Returns the list of features within the specified reference sequence and range """
        args = range_arguments.parse_args(request)
        start = args.get('start')
        end = args.get('end')

        fq = db.session.query(Feature).join(RefSeq).join(Species).filter(RefSeq.locus == ref_seq, Species.name == species, Feature.end >= start, Feature.start <= end)

        features = []

        for gene in fq.filter(Feature.feature == 'gene').all():
            gene_feature = enumerate_feature(gene)
            self.add_subfeatures(gene_feature, fq)
            features.append(gene_feature)

        return {'features': features}


@ns.route('/jbrowse/<string:species>/stats/global')
class JbrowseGLobalStatsAPI(Resource):
    def get(self, species):
        """ Returns global statistics about features served from this store """

        n_features = 0
        seq_bps = 0

        for ref_seq in db.session.query(RefSeq).join(Species).filter(Species.name == species).all():
            n_features = db.session.query(Feature).join(RefSeq).join(Species).filter(Species.name == species, RefSeq.name == ref_seq.name).count()
            seq_bps += ref_seq.length

        return {'featureDensity': float(n_features)/seq_bps, 'featureCount': n_features}


@ns.route('/jbrowse/<string:species>/stats/region/<string:ref_seq>')
@api.response(404, 'Reference sequence not found.')
class JbrowseRegionalStatsAPI(Resource):
    @api.expect(range_arguments, validate=True)
    def get(self, species, ref_seq):
        """ Returns statistics for a particular region """
        args = range_arguments.parse_args(request)
        start = args.get('start')
        end = args.get('end')

        seq = db.session.query(RefSeq).join(Species).filter(Species.name == species, RefSeq.name == ref_seq).one_or_none()

        if not seq:
            return None, 404

        n_features = db.session.query(Feature).join(RefSeq).join(Species).filter(Species.name == species, RefSeq.name == seq.name, Feature.end >= start, Feature.start <= end).count()
        return {'featureDensity': float(n_features)/seq.length, 'featureCount': n_features}


@ns.route('/sequences/<string:species>/<string:ref_seq>')
@api.response(404, 'Reference sequence not found.')
class SequencesAPI(Resource):
    def get(self, species, ref_seq):
        """ Returns nucleotide sequences from selected references. Use 'all' to wildcard species or ref_seq """

        ref_seqs = db.session.query(RefSeq)

        if species != 'all':
            ref_seqs = ref_seqs.filter(Species.name == species)

        if ref_seq != 'all':
            ref_seqs = ref_seqs.filter(RefSeq.name == ref_seq)

        if not ref_seqs.count():
            return None, 404

        sequences = []

        for ref in ref_seqs.all():
            x_keys = {}
            if species == 'all':
                x_keys['species'] = ref.species.name
            if ref_seq == 'all':
                x_keys['ref_seq'] = ref.name
            for f in ref.features:
                if f.sequence:
                    seq = object_as_dict(f.sequence)
                    del seq['id']
                    del seq['species_id']
                    if len(x_keys):
                        seq.update(x_keys)
                    sequences.append(seq)

        return {'sequences': sequences}
