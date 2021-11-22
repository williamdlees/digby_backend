from flask import Flask, render_template, request, flash, Blueprint, redirect, url_for
from flask_migrate import Migrate
from flask_security import Security, SQLAlchemyUserDatastore, login_required
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_admin import Admin
from flask_cors import CORS
import os
import custom_logging
from reverse_proxied import ReverseProxied
from db.vdjbase_db import vdjbase_db_init, manage_airrseq, airrseq_import, airrseq_copy, airrseq_remove

from flask_security.utils import hash_password
from flask_sqlalchemy import SQLAlchemy

from extensions import celery

sql_db = None


app = Flask(__name__)
bootstrap = Bootstrap(app)
app.config.from_pyfile('config.cfg')
app.config.from_pyfile('secret.cfg')

# configure/initialize all your extensions

sql_db = celery.init_app(app)



app.wsgi_app = ReverseProxied(app.wsgi_app)

CORS(app)


app.config['BASE_PATH'] = os.getcwd()

if 'STATIC_PATH' not in app.config:
    app.config['STATIC_PATH'] = os.path.join(app.config['BASE_PATH'], 'static')

if 'OUTPUT_PATH' not in app.config:
    app.config['OUTPUT_PATH'] = os.path.join(app.config['STATIC_PATH'], 'output')

if 'UPLOAD_PATH' not in app.config:
    app.config['UPLOAD_PATH'] = os.path.join(app.config['BASE_PATH'], 'uploads')

app.config['R_SCRIPT_PATH'] = os.path.join(app.config['BASE_PATH'], 'api/reports/R_scripts')

if 'R_LIBS' not in os.environ or os.environ['R_LIBS'] is None or len(os.environ['R_LIBS']) < 1:
    os.environ['R_LIBS'] = app.config['R_SCRIPT_PATH']

# TODO - make this work on Windows as well as Unix, the D:\ screws it up
"""
else:
    r_libs = os.environ['R_LIBS'].split(':')
    r_libs = r_libs.append(app.config['R_SCRIPT_PATH'])
    os.environ['R_LIBS'] = ':'.join(r_libs)
"""

custom_logging.init_logging(app)

mail = Mail(app)


vdjbase_dbs = vdjbase_db_init(os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db'))

admin_obj = Admin(app, template_mode='bootstrap3')

from security.useradmin import *
from security.security import *

user_datastore = SQLAlchemyUserDatastore(db, User, Role)
security = Security(app, user_datastore, confirm_register_form=ExtendedRegisterForm)

from api.restx import api
from api.genomic.genomic import ns as genomic
from api.vdjbase.vdjbase import ns as vdjbase
from api.reports.reports import ns as reports

from db.feature_db import *
from db.update import update_genomic_db
from db.build_gff import build_gffs
import db.vdjbase_maint
import db.vdjbase_export
from db.vdjbase_igsnper import do_igsnper

migrate = Migrate(app, sql_db)

blueprint = Blueprint('api', __name__, url_prefix='/api')
api.init_app(blueprint)
api.add_namespace(genomic)
api.add_namespace(vdjbase)
api.add_namespace(reports)
app.register_blueprint(blueprint)

from api.reports.reports import load_report_defs
load_report_defs()

app.logger.error('Digby backend started')


@app.route('/', methods=['GET', 'POST'])
def index():
    if user_datastore.find_role('Admin') is None:
        return redirect(url_for('create_user'))

    return render_template('index.html', current_user=current_user)


@app.route('/create_user', methods=['GET', 'POST'])
def create_user():
    if user_datastore.find_role('Admin') is not None:
        return redirect('/')

    form = FirstAccountForm()

    if request.method == 'POST':
        if form.validate():
            user = user_datastore.create_user(email=form.email.data, password=hash_password(form.password.data), name=form.name.data)
            sql_db.session.commit()
            user_datastore.create_role(name='Admin')
            user_datastore.add_role_to_user(user, 'Admin')
            sql_db.session.commit()
            flash("User created")
            return redirect('/')

    return render_template('security/first_account.html', form=form)


@app.route('/profile', methods=['GET', 'POST'])
@login_required
def profile():
    form = ProfileForm(obj=current_user)
    form.email = ''
    if request.method == 'POST':
        if form.validate():
            save_Profile(db, current_user, form)
            flash('Profile updated.')

    return render_template('profile.html', form=form, current_user=current_user, url='profile')


@app.route('/update_genomic', methods=['GET', 'POST'])
@login_required
def update_genomic():
    return update_genomic_db()


@app.route('/airrseq', methods=['GET', 'POST'])
def airrseq():
    return manage_airrseq(app)


@app.route('/airrseq_import_status/<species>/<dataset>', methods=['GET', 'POST'])
def airrseq_import_status(species, dataset):
    return airrseq_import(species.replace(' ', '_'), dataset.replace(' ', '_'), app)


@app.route('/airrseq_copy_live', methods=['GET', 'POST'])
def airrseq_copy_live():
    return airrseq_copy(app, vdjbase_dbs)


@app.route('/airrseq_delete/<species>/<dataset>', methods=['GET', 'POST'])
def airrseq_delete(species, dataset):
    return airrseq_remove(species, dataset, app, vdjbase_dbs)


@app.route('/build_gff', methods=['GET', 'POST'])
@login_required
def build_gff():
    return build_gffs()


@app.route('/export_vdjbase_metadata', methods=['GET', 'POST'])
@login_required
def export_vdjbase_metadata():
    return db.vdjbase_export.export_metadata()


@app.route('/create_igsnp/<species>/<dataset>/', methods=['GET', 'POST'])
@login_required
def create_igsnp(species, dataset):
    return do_igsnper(species, dataset)
