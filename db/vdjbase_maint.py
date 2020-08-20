#
# backend maintenance functions for VDJbase-style databases
#
from app import app, vdjbase_dbs
import os.path
import os
import shutil

from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool

from db.vdjbase_confidence import check_novel_confidence
from db.vdjbase_cross_sample import update_alleles_appearance, calculate_gene_frequencies, calculate_patterns
from db.vdjbase_model import *
from db.vdjbase_db import Session
from db.vdjbase_reference import import_reference_alleles
from db.vdjbase_projects import import_studies, process_genotypes, add_deleted_alleles, process_haplotypes_and_stats
from db.vdjbase_exceptions import *


#
# Create a sqlite database using the files in the export directory's vdjbase_metadata subdirectory. The
# database is put in that subdirectory and can then be moved by hand to the appropriate location in static
#
def create_databases():
    export_dir = os.path.join(app.config['EXPORT_DIR'], 'vdjbase_metadata')
    results = []

    if not os.path.isdir(export_dir):
        return 'vdjbase_metadata folder does not exist - aborting'

    for s_item in os.listdir(export_dir):
        if os.path.isdir(os.path.join(export_dir, s_item)) and s_item[0] != '.':
            for d_item in os.listdir(os.path.join(export_dir, s_item)):
                if os.path.isdir(os.path.join(export_dir, s_item, d_item)) and d_item[0] != '.':
                    results.append(create_single_database(export_dir, s_item, d_item))

    return '\n'.join(results)

def create_single_database(export_dir, species, dataset):
    result = ['Processing %s/%s' %(species, dataset)]
    ds_dir = os.path.join(export_dir, species, dataset)

    if not os.path.isdir(ds_dir):
        result.append('data folder for %s/%s does not exist - skipped' % (species, dataset))
        return '\r\n'.join(result)

    reference_dir = os.path.join(ds_dir, 'reference')
    if not os.path.isdir(reference_dir):
        result.append('Reference directory for %s/%s not found - skipped' % (species, dataset))
        return '\r\n'.join(result)

    db_file = os.path.join(ds_dir, 'db.sqlite3')

    if os.path.isfile(db_file):
        os.remove(db_file)
    # shutil.copyfile(os.path.join(ds_dir, 'copy_db.sqlite3'), db_file)

    engine = create_engine('sqlite:///' + db_file, echo=False, poolclass=NullPool)
    Base.metadata.create_all(engine)
    db_connection = engine.connect()
    engine.session = Session(bind=db_connection)
    session = engine.session

    try:
        result.extend(import_reference_alleles(reference_dir, session, species))
        result.extend(import_studies(ds_dir, species, dataset, session))
        result.extend(add_deleted_alleles(session))
        result.extend(process_genotypes(ds_dir, species, dataset, session))
        result.extend(process_haplotypes_and_stats(ds_dir, species, dataset, session))
        result.extend(update_alleles_appearance(session))
        result.extend(calculate_gene_frequencies(ds_dir, session))
        result.extend(calculate_patterns(session))
        result.extend(check_novel_confidence(ds_dir, session))
    except DbCreationError as e:
        result = [e.args[0]]
    finally:
        db_connection.close()
        engine.dispose()

    return '<br>'.join(result)





