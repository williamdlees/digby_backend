# Standalone script to make a genomic sqlite database from the command line

import argparse
from db.genomic_maint import create_dataset
import os

parser = argparse.ArgumentParser(description='Make a genomic sqlite database from files in current directory')
args = parser.parse_args()

create_dataset(os.getcwd())