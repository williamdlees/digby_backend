# Standalone script to maake a VDJbase sqlite database from the command line

import argparse
from db.vdjbase_maint import create_single_database
import os

parser = argparse.ArgumentParser(description='Make a VDJbase sqlite database from files n current directory')
parser.add_argument('species', help='species')
parser.add_argument('dataset_name', help='data set name')
args = parser.parse_args()

class Job:
    def update_state(self, meta="", state=""):
        print('meta: %s state: %s' % (meta['value'], state))


create_single_database(Job(), args.species, args.dataset_name, os.getcwd(), True)
