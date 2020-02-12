# Long read sequences submitted to IMGT

from flask_restplus import Resource
from api.restplus import api
import csv
import json

ns = api.namespace('long_read/submitted_sequences', description='Long read sequences submitted to IMGT')


@ns.route('/')
class SubmittedSequences(Resource):

    def get(self):
        """
        Returns list of sequences submitted to IMGT
        """
        data = []

        for row in csv.DictReader(open('db/submitted.csv', 'r', encoding='ISO-8859-1')):
            data.append(row)

        return json.dumps(data)
