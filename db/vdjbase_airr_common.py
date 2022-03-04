# Common utilities for managing miairr integration

import csv
import os

definition_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vdjbase_airr_schema_defs.csv')


def read_definition_data(remove_dots=True):
    defs = {}
    with open(definition_file, 'r') as fi:
        reader = csv.DictReader(fi)
        for row in reader:
            if 'VDJbase table' in row and row['VDJbase table']:
                if row['VDJbase table'] not in defs:
                    defs[row['VDJbase table']] = {}

                if remove_dots:
                    row['simple_name'] = row['simple_name'].replace('.', '_')

                defs[row['VDJbase table']][row['simple_name']] = row

    return defs


