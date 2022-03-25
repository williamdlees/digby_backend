# Standalone script to make a genomic sqlite database from the command line

import argparse

from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from sqlalchemy.pool import NullPool

from db.genomic_db import Base
from db.genomic_maint import create_dataset
from db.build_gff import build_gff
import os

parser = argparse.ArgumentParser(description='Make a genomic sqlite database from files in current directory')
parser.add_argument('species', help='species')
parser.add_argument('dataset_name', help='data set name')
args = parser.parse_args()

#create_dataset(args.species, args.dataset_name)
#quit()

# or comment out the above lines to build the gffs without rebuilding the database

dataset_dir = os.getcwd()
db_file = os.path.join(dataset_dir, 'db.sqlite3')
engine = create_engine('sqlite:///' + db_file, echo=False, poolclass=NullPool)
Base.metadata.create_all(engine)
db_connection = engine.connect()
engine.session = Session(bind=db_connection)
session = engine.session

build_gff(session, dataset_dir)
