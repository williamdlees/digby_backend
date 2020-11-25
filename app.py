from flask import Flask, render_template, request, flash, Blueprint
from flask_migrate import Migrate
from flask_security import Security, SQLAlchemyUserDatastore, login_required
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_admin import Admin
from flask_cors import CORS
import os
import logging.handlers
from reverse_proxied import ReverseProxied

from flask_security.utils import hash_password
from flask_sqlalchemy import SQLAlchemy

from extensions import celery

sql_db = None


def create_app():
    global sql_db
    app = Flask(__name__)
    app.config.from_pyfile('config.cfg')
    app.config.from_pyfile('secret.cfg')

    # configure/initialize all your extensions

    sql_db = celery.init_app(app)

    return app


app = create_app()
app.wsgi_app = ReverseProxied(app.wsgi_app)

CORS(app)

bootstrap = Bootstrap(app)

app.config['BASE_PATH'] = os.getcwd()

if 'STATIC_PATH' not in app.config:
    app.config['STATIC_PATH'] = os.path.join(app.config['BASE_PATH'], 'static')

if 'OUTPUT_PATH' not in app.config:
    app.config['OUTPUT_PATH'] = os.path.join(app.config['STATIC_PATH'], 'output')

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

handler = logging.handlers.RotatingFileHandler('app.log', maxBytes=1024 * 1024)
handler.setLevel(logging.INFO)
app.logger.addHandler(handler)

mail = Mail(app)

from db.vdjbase_db import vdjbase_db_init

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


@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', current_user=current_user)


@app.route('/create_user', methods=['GET', 'POST'])
def create_user():
    if user_datastore.find_role('Admin') is None:                       # First live user gets admin rights
        user = user_datastore.create_user(email="admin@vdjbase.org", password=hash_password("admin"))
        sql_db.session.commit()
        user_datastore.create_role(name='Admin')                        # You will want to remove this in any real application!
        user_datastore.add_role_to_user(user, 'Admin')          # You will want to remove this in any real application!
        sql_db.session.commit()
        return "User created"
    else:
        return "User not created"


@app.route('/profile', methods=['GET', 'POST'])
@login_required
def profile():
    form = ProfileForm(obj = current_user)
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


@app.route('/build_gff', methods=['GET', 'POST'])
@login_required
def build_gff():
    return build_gffs()


@app.route('/export_vdjbase_metadata', methods=['GET', 'POST'])
@login_required
def export_vdjbase_metadata():
    return db.vdjbase_export.export_metadata()


@app.route('/create_vdjbase_db', methods=['GET', 'POST'])
@login_required
def create_vdjbase_db():
    return db.vdjbase_maint.create_databases()


@app.route('/create_igsnp/<species>/<dataset>/', methods=['GET', 'POST'])
@login_required
def create_igsnp(species, dataset):
    return do_igsnper(species, dataset)
