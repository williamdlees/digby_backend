# IMGT Gene Table

from flask_restplus import Resource
from api.restplus import api
import json
from api.fetch_imgt_alleles import fetch_imgt_alleles

ns = api.namespace('imgt/imgt_gene_table', description='IMGT Gene Table')


@ns.route('/')
class ImgtAlleles(Resource):

    def get(self):
        """
        Returns IMGT Gene Table
        """
        data = fetch_imgt_alleles('human', 'IGHV')
        return json.dumps(data)
