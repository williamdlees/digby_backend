# Manage a list of available vdjbase-style databases

from os.path import join, isdir
from os import listdir
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from app import app
import os.path

VDJBASE_DB_PATH = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db')

Session = sessionmaker()

class ContentProvider():
    db = None
    connection = None
    session = None

    def __init__(self, path):
        self.db = create_engine('sqlite:///' + path + '?check_same_thread=false', echo=False, pool_threadlocal=True)
        self.connection = self.db.connect()
        self.session = Session(bind=self.connection)


def vdjbase_db_init():
    vdjbase_dbs = {}

    for species in listdir(VDJBASE_DB_PATH):
        p = join(VDJBASE_DB_PATH, species)
        if isdir(p) and species[0] != '.':
            for name in listdir(p):
                if isdir(join(p, name)) and name[0] != '.':
                    if species not in vdjbase_dbs:
                        vdjbase_dbs[species] = {}
                    vdjbase_dbs[species][name] = ContentProvider(join(p, name, 'db.sqlite3'))

    return vdjbase_dbs

