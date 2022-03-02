# Common utilities for managing miairr integration

import csv

definition_file = 'vdjbase_airr_schema_defs.csv'


def read_definition_data():
    defs = {}
    with open(definition_file, 'r') as fi:
        reader = csv.DictReader(fi)
        for row in reader:
            if 'VDJbase table' in row and row['VDJbase table']:
                if row['VDJbase table'] not in defs:
                    defs[row['VDJbase table']] = {}

                row['simple_name'] = row['simple_name'].replace('.', '_')
                defs[row['VDJbase table']][row['simple_name']] = row

    return defs


