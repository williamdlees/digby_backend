import os
import custom_logging
import yaml

from werkzeug.middleware.proxy_fix import ProxyFix
from flask import Flask, render_template, Blueprint, redirect, url_for, send_from_directory
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_admin import Admin
from flask_cors import CORS
from flask_swagger_ui import get_swaggerui_blueprint
from extensions import celery

from db.vdjbase_db import study_data_db_init
from db.madc import madc_init

sql_db = None

app = Flask(__name__, static_folder=None)
app.static_url_path = '/static'
bootstrap = Bootstrap(app)
app.static_url_path = None
app.config.from_pyfile('config.cfg')
app.config.from_pyfile('secret.cfg')

# configure/initialize all your extensions

sql_db = celery.init_app(app)

app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)

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


mail = Mail(app)
custom_logging.init_logging(app, mail)

vdjbase_dbs = study_data_db_init(os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db'))
genomic_dbs = study_data_db_init(os.path.join(app.config['STATIC_PATH'], 'study_data/Genomic/db'))
madc_index = madc_init(app)

admin_obj = Admin(app, template_mode='bootstrap3')

from api.restx import api
from api.genomic.genomic import ns as genomic
from api.vdjbase.vdjbase import ns as vdjbase
from api.reports.reports import ns as reports
from api.system.system import ns as system, digby_protected

from db.genomic_db import *
import db.vdjbase_maint
import db.vdjbase_export
from db.vdjbase_igsnper import do_igsnper


blueprint = Blueprint('api', __name__, url_prefix='/api')
api.init_app(blueprint)
api.add_namespace(genomic)
api.add_namespace(vdjbase)
api.add_namespace(reports)
api.add_namespace(system)
app.register_blueprint(blueprint)

from api_v1.open_api import api_bp

# PROXY_PREFIX is needed to get around a bug in swagger_ui: https://github.com/PWZER/swagger-ui-py/issues/34
# When this is fixed, ProxyFix will take care of things and the special location for /admin/api_v1 in nginix confs can be removed

app.register_blueprint(api_bp, url_prefix="/api/v1")

SWAGGER_URL = "/api/v1"
API_URL = '/static/vdjbase-api-openapi3.yaml'
swaggerui_blueprint = get_swaggerui_blueprint(
    SWAGGER_URL,
    API_URL,
    config={
        'app_name': "VDJbase Open API v1"
    }
)


@app.route('/schema/<path:filename>')
def serve_schema(filename):
    return send_from_directory('schema', filename)


app.register_blueprint(swaggerui_blueprint, url_prefix=SWAGGER_URL)

with open('schema/vdjbase-api-openapi3.yaml', 'r') as f:
    openapi_schema = yaml.safe_load(f)

from api.reports.reports import load_report_defs
load_report_defs()

from flask_jwt_extended import JWTManager

app.config["JWT_TOKEN_LOCATION"] = ["headers"]
jwt = JWTManager(app)


@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', current_user='')


@app.route('/static/<path:path>', methods=['GET', 'POST'])
def static(path):
    if '/gff' in path:
        return send_from_gff(path)
    else:
        return send_from_directory(app.config['STATIC_PATH'], path)

@digby_protected()
def send_from_gff(path):
    return send_from_directory(app.config['STATIC_PATH'], path)




#@app.route('/export_vdjbase_metadata', methods=['GET', 'POST'])
#@login_required
#def export_vdjbase_metadata():
#    return db.vdjbase_export.export_metadata()


#@app.route('/create_igsnp/<species>/<dataset>/', methods=['GET', 'POST'])
#@login_required
#def create_igsnp(species, dataset):
#    return do_igsnper(species, dataset)
