# Manage a list of available vdjbase-style databases

from os.path import join, isdir, isfile
from os import listdir
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from app import app
import os.path


Session = sessionmaker()

class ContentProvider():
    db = None
    connection = None
    session = None

    def __init__(self, path):
        self.db = create_engine('sqlite:///' + path + '?check_same_thread=false', echo=False, pool_threadlocal=True)
        self.connection = self.db.connect()
        self.session = Session(bind=self.connection)

    def close(self):
        self.connection.close()
        self.db.dispose()


def vdjbase_db_init(vdjbase_db_path):
    vdjbase_dbs = {}

    for species in listdir(vdjbase_db_path):
        p = join(vdjbase_db_path, species)
        if isdir(p) and species[0] != '.':
            for name in listdir(p):
                if isdir(join(p, name)) and name[0] != '.' and '.txt' not in name:
                    description = ''
                    if isfile(join(p, name, 'db_description.txt')):
                        with open(join(p, name, 'db_description.txt'), 'r') as fi:
                            description = ' '.join(fi.readlines())
                    if species not in vdjbase_dbs:
                        vdjbase_dbs[species] = {}
                    vdjbase_dbs[species][name] = ContentProvider(join(p, name, 'db.sqlite3'))
                    vdjbase_dbs[species][name + '_description'] = description
    return vdjbase_dbs


