# Manage a list of available vdjbase-style databases
import os
import shutil
from os.path import join, isdir, isfile
from os import listdir
from time import sleep
from flask import render_template, request, redirect, url_for, Markup
from sqlalchemy import create_engine, inspect, text
from sqlalchemy.orm import sessionmaker
from flask_table import Table, Col
from flask_wtf import FlaskForm
from werkzeug.exceptions import BadRequest
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from flask_wtf.file import FileField, FileRequired
from wtforms.validators import DataRequired
from werkzeug.utils import secure_filename

from db.vdjbase_exceptions import DbCreationError
from db.vdjbase_maint import create_single_database
from extensions import celery
import traceback

Session = sessionmaker()

class ContentProvider():
    db = None
    connection = None
    session = None

    def __init__(self, path):
        self.db = create_engine('sqlite:///' + path + '?check_same_thread=false', echo=False)
        self.connection = self.db.connect()
        self.session = Session(bind=self.connection)

    def close(self):
        self.connection.close()
        self.db.dispose()


# TODO: when we change to biomial species names, update this lookup and corresponding code
species_lookup = {
    'Human': ('Homo sapiens', 'NCBITAXON: 9606'),
    'Rhesus Macaque': ('Macaca mulatta', 'NCBITAXON: 9544'),
    'Crab-eating Macaque': ('Macaca fascicularis', 'NCBITAXON: 9541'),
}


def study_data_db_init(vdjbase_db_path):
    sqlite_dbs = {}

    for species in listdir(vdjbase_db_path):
        p = join(vdjbase_db_path, species)
        if isdir(p) and species[0] != '.':
            for name in listdir(p):
                if isdir(join(p, name)) and name[0] != '.' and '.txt' not in name:
                    description = ''
                    if isfile(join(p, name, 'db_description.txt')):
                        with open(join(p, name, 'db_description.txt'), 'r') as fi:
                            description = ' '.join(fi.readlines())
                    if species not in sqlite_dbs:
                        sqlite_dbs[species] = {}
                    sqlite_dbs[species][name] = ContentProvider(join(p, name, 'db.sqlite3'))
                    sqlite_dbs[species][name + '_description'] = {'description_text': description}

                    if species in species_lookup:
                        sqlite_dbs[species][name + '_description']['binomial'] = species_lookup[species][0]
                        sqlite_dbs[species][name + '_description']['taxid'] = species_lookup[species][1]
                    else:
                        sqlite_dbs[species][name + '_description']['binomial'] = species
                        sqlite_dbs[species][name + '_description']['taxid'] = ''
                        print('Species %s not found in species lookup' % species)

    # temp fix: add asc_genotype column to sample table if not there already

    if 'genomic' not in vdjbase_db_path.lower():
        for species in sqlite_dbs:
            for locus in sqlite_dbs[species]:
                if 'description' not in locus:
                    inspector = inspect(sqlite_dbs[species][locus].db)
                    cols = inspector.get_columns('Sample')
                    if 'asc_genotype' not in [col['name'] for col in cols]:
                        with sqlite_dbs[species][locus].connection as con:
                            con.execute('ALTER TABLE Sample ADD COLUMN asc_genotype text')
                            sqlite_dbs[species][locus].session.commit()


    # sort datasets of each species

    for species in sqlite_dbs:
        sqlite_dbs[species] = dict(sorted(sqlite_dbs[species].items(), key=lambda kv: kv[0]))

    # put Human at the front
    sqlite_dbs = dict(sorted(sqlite_dbs.items(), key=lambda kv: 'aaaaa' if kv[0] == 'Human' else kv[0]))

    return sqlite_dbs




def manage_airrseq(app):
    class ItemTable(Table):
        species = Col('Species')
        dataset = Col('Dataset')
        action = Col('Action')
        classes = ['table']

    class Item(object):
        def __init__(self, species, dataset):
            self.species = species
            self.dataset = dataset
            self.action = Markup('<a class="btn" href="airrseq_delete/%s/%s">Delete</a>') % (species, dataset)

    class AddDatasetForm(FlaskForm):
        species = StringField('Species', validators=[DataRequired()])
        dataset = StringField('Dataset', validators=[DataRequired()])
        datafile = FileField(validators=[FileRequired()])
        submit = SubmitField('Add New Dataset')

    vdjbase_db_path = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db')
    items = []

    for species in listdir(vdjbase_db_path):
        p = join(vdjbase_db_path, species)
        if isdir(p) and species[0] != '.':
            for name in listdir(p):
                if isdir(join(p, name)) and name[0] != '.' and '.txt' not in name:
                    items.append(Item(species, name))

    form = AddDatasetForm()

    if form.validate_on_submit():
        upload_path = app.config['UPLOAD_PATH']
        species = secure_filename(form.species.data).replace(' ', '_')
        dataset = secure_filename(form.dataset.data).replace(' ', '_')

        if not os.path.exists(os.path.join(upload_path, species)):
            os.mkdir(os.path.join(upload_path, species))

        dataset_path = os.path.join(upload_path, species, dataset)

        if os.path.exists(dataset_path):
            replacing = 'Replacing existing dataset %s - %s' % (species, dataset)
            shutil.rmtree(dataset_path, ignore_errors=True)

        os.mkdir(dataset_path)
        data = request.files[form.datafile.name].read()
        open(os.path.join(dataset_path, 'data.zip'), 'wb').write(data)
        return redirect(url_for('airrseq_import_status', species=species, dataset=dataset))

    return render_template('airrseq.html', item_table=ItemTable(items), form=form)


class Job:
    def update_state(self, meta="", state=""):
        print('meta: %s state: %s' % (meta['value'], state))


def airrseq_import(species, dataset, app):
    # result = airrseq_import_celery.delay(species, dataset, app.config['UPLOAD_PATH'], app.config['STATIC_PATH'])
    create_single_database(Job(), species, dataset, app.config['UPLOAD_PATH'])
    return ""
    return render_template('airrseq_import_status.html', species=species, dataset=dataset, id=result.id)


@celery.task(bind=True)
def airrseq_import_celery(self, species, dataset, upload_path, static_path):
    try:
        return import_airrseq_data(self, species, dataset, upload_path, static_path)
    except BadRequest as bad:
        print('BadRequest raised during report processing: %s' % bad.description)
        return {'status': 'error', 'description': bad.description}
    except Exception as e:
        print('Exception raised during report processing: %s' % traceback.format_exc())
        return {'status': 'error', 'description': 'Unexpected error when importing airrseq dataset: %s' % traceback.format_exc()}


def import_airrseq_data(job, species, dataset, upload_path, static_path):
    status, results = create_single_database(job, species, dataset, upload_path)

    if status:
        results.insert(0, 'Dataset created!<br><br>')
    else:
        results.insert(0, 'Dataset not created - see logs below:<br><br>')

    return {'status': 'ok', 'log': '<br>'.join(results)}


def airrseq_copy(app, vdjbase_dbs):
    class CopyDatasetForm(FlaskForm):
        species = StringField('Species', validators=[DataRequired()])
        dataset = StringField('Dataset', validators=[DataRequired()])
        description = StringField('Dataset Description', validators=[DataRequired()])
        submit = SubmitField('Make Dataset Live')

    form = CopyDatasetForm()

    if form.validate_on_submit():
        species = secure_filename(form.species.data).replace(' ', '_')
        dataset = secure_filename(form.dataset.data).replace(' ', '_')

        try:
            do_airrseq_copy(species, dataset, form.description.data, app, vdjbase_dbs)
        except DbCreationError as e:
            print('Exception encountered processing report request: %s' % traceback.format_exc())
            app.logger.error('Exception encountered processing report request: %s' % traceback.format_exc())
            return render_template('airrseq_copy.html', form=form, error='Exception encountered processing report request: %s' % traceback.format_exc())

        return redirect(url_for('airrseq'))

    return render_template('airrseq_copy.html', form=form, error=None)


def do_airrseq_copy(species, dataset, description, app, vdjbase_dbs):
    our_upload_path = os.path.join(app.config['UPLOAD_PATH'], species, dataset)

    if not os.path.isdir(our_upload_path) or not os.path.isfile(os.path.join(our_upload_path, 'db.sqlite3')):
        raise DbCreationError('Dataset %s/%s not found' % (species, dataset))

    if not os.path.isdir(os.path.join(our_upload_path, 'samples')):
        raise DbCreationError('Sample directory for %s/%s not found' % (species, dataset))

    try:
        if species in vdjbase_dbs and dataset in vdjbase_dbs[species]:
            vdjbase_dbs[species][dataset].close()
            del(vdjbase_dbs[species][dataset])

        vdjbase_db_path = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db')
        vdjbase_sample_path = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/samples')
        our_db_path = os.path.join(vdjbase_db_path, species, dataset)
        our_sample_path = os.path.join(vdjbase_sample_path, species, dataset)

        if not os.path.isdir(os.path.join(vdjbase_db_path, species)):
            os.mkdir(os.path.join(vdjbase_db_path, species))

        if not os.path.isdir(os.path.join(vdjbase_sample_path, species)):
            os.mkdir(os.path.join(vdjbase_db_path, species))


        if isdir(our_db_path):
            shutil.rmtree(our_db_path, ignore_errors=True)

        if isdir(our_sample_path):
            shutil.rmtree(our_sample_path, ignore_errors=True)

        os.mkdir(our_db_path)
        os.mkdir(our_sample_path)

        shutil.copyfile(os.path.join(our_upload_path, 'db.sqlite3'), os.path.join(our_db_path, 'db.sqlite3'))
        with open(os.path.join(our_db_path, 'db_description.txt'), 'w') as fo:
            fo.write(description)

        for node in os.listdir(os.path.join(our_upload_path, 'samples')):
            if node[0] != '.' and os.path.isdir(os.path.join(our_upload_path, 'samples', node)):
                shutil.copytree(os.path.join(our_upload_path, 'samples', node), os.path.join(our_sample_path, node))

        if species not in vdjbase_dbs:
            vdjbase_dbs[species] = {}

        vdjbase_dbs[species][dataset] = ContentProvider(os.path.join(our_db_path, 'db.sqlite3'))
        vdjbase_dbs[species][dataset + '_description'] = description

    except Exception as e:
        raise DbCreationError(str(e))


def airrseq_remove(species, dataset, app, vdjbase_dbs):
    vdjbase_db_path = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db')
    vdjbase_sample_path = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/samples')
    our_db_path = os.path.join(vdjbase_db_path, species, dataset)
    our_sample_path = os.path.join(vdjbase_sample_path, species, dataset)

    if species in vdjbase_dbs and dataset in vdjbase_dbs[species]:
        vdjbase_dbs[species][dataset].close()
        del(vdjbase_dbs[species][dataset])

    if isdir(our_db_path):
        shutil.rmtree(our_db_path, ignore_errors=True)

    if isdir(our_sample_path):
        shutil.rmtree(our_sample_path, ignore_errors=True)

    return redirect(url_for('airrseq'))
