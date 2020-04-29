from flask import Flask, render_template, request, flash, Blueprint
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_security import Security, SQLAlchemyUserDatastore, login_required
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_admin import Admin
from flask_cors import CORS

import logging.handlers

print("app name: %s" % __name__)
app = Flask(__name__)
CORS(app)

bootstrap = Bootstrap(app)
app.config.from_pyfile('config.cfg')
app.config.from_pyfile('secret.cfg')

handler = logging.handlers.RotatingFileHandler('app.log', maxBytes=1024 * 1024)
handler.setLevel(logging.INFO)
app.logger.addHandler(handler)

mail = Mail(app)

db = SQLAlchemy(app)

admin_obj = Admin(app, template_mode='bootstrap3')

from security.useradmin import *
from security.security import *

user_datastore = SQLAlchemyUserDatastore(db, User, Role)
security = Security(app, user_datastore, confirm_register_form=ExtendedRegisterForm)

from api.restplus import api
from api.genomic.genomic import ns as genomic

from db.feature_db import *
from db.update import db_update
from db.build_gff import build_gffs

migrate = Migrate(app, db)

blueprint = Blueprint('api', __name__, url_prefix='/api')
api.init_app(blueprint)
api.add_namespace(genomic)
app.register_blueprint(blueprint)

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', current_user=current_user)


@app.route('/profile', methods=['GET', 'POST'])
@login_required
def profile():
    if user_datastore.find_role('Admin') is None:                       # First live user gets admin rights
        user_datastore.create_role(name='Admin')                        # You will want to remove this in any real application!
        user_datastore.add_role_to_user(current_user, 'Admin')          # You will want to remove this in any real application!
        db.session.commit()                                             # You will want to remove this in any real application!

    form = ProfileForm(obj = current_user)
    form.email = ''
    if request.method == 'POST':
        if form.validate():
            save_Profile(db, current_user, form)
            flash('Profile updated.')

    return render_template('profile.html', form=form, current_user=current_user, url='profile')


@app.route('/update', methods=['GET', 'POST'])
@login_required
def update():
    return db_update()

@app.route('/build', methods=['GET', 'POST'])
@login_required
def build():
    return build_gffs()


@app.route('/play', methods=['GET', 'POST'])
@login_required
def play():
    feats = db.session.query(Feature).all()

    for feat in feats:
        for sam in feat.samples:
            print(sam.sample.type)
