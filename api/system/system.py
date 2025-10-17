# Services related to vdjbase repseq-based data sets
from functools import wraps

import requests
from flask import jsonify
from flask_jwt_extended import create_access_token, set_access_cookies, jwt_required, get_jwt_identity, \
    verify_jwt_in_request, decode_token, create_refresh_token
from flask_restx import Resource, reqparse, fields, marshal, inputs
from api.restx import api
import json
from app import vdjbase_dbs, app
from datetime import datetime
import time


ns = api.namespace('system', description='System-specific information')


@ns.route('/config')
class ConfigApi(Resource):
    def get(self):
        """ Return internal configuration details """

        config = {
            'app_protected': not (app.config['JWT_USER'] == '' and app.config['JWT_PASSWORD'] == ''),
            'wp_news_url': app.config['WORDPRESS_NEWS_URL'],
            'wp_help_url': app.config['WORDPRESS_HELP_URL'],
            'wp_rest': app.config['WORDPRESS_REST']}

        # get category details from Wordpress
        #    https://renemorozowich.com/using-wordpress-rest-api-get-blogs/

        wanted_slugs = ['vdjbase_news', 'vdjbase_help']
        wp_url = {}
        wp_url['vdjbase_news'] = app.config['WORDPRESS_NEWS_URL'] + app.config['WORDPRESS_REST']
        wp_url['vdjbase_help'] = app.config['WORDPRESS_HELP_URL'] + app.config['WORDPRESS_REST']

        for wanted_slug in wanted_slugs:
            r = requests.get(wp_url[wanted_slug] + 'categories')
            if r.status_code == 200:
                resp = r.content.decode("utf-8")
                resp = json.loads(resp)

                for rec in resp:
                    if rec['slug'] == wanted_slug:
                        config[rec['slug']] = '%sposts?categories=%s' % (wp_url[wanted_slug], rec['id'])

        return config


@ns.route("/login/<username>/<password>")
class LoginApi(Resource):
    def get(self, username, password):
        if username == app.config['JWT_USER'] and password == app.config['JWT_PASSWORD']:
            access_token = create_access_token(identity=app.config['JWT_USER'])
            access_token_details = decode_token(access_token)
            refresh_token = create_refresh_token(identity=app.config['JWT_USER'])
            refresh_token_details = decode_token(refresh_token)
            response = jsonify({
                "msg": "login successful",
                "username": username,
                "access_token": access_token,
                "access_token_lifetime": access_token_details['exp'] - int(time.time()),
                "refresh_token": refresh_token,
                "refresh_token_lifetime": refresh_token_details['exp'] - int(time.time()),
            })
            set_access_cookies(response, access_token)
            return response

        return "Unrecognised user name or password", 401



@ns.route("/refresh")
class RefreshApi(Resource):
    @jwt_required(refresh=True)
    def get(self):
        identity = get_jwt_identity()
        access_token = create_access_token(identity=identity)
        access_token_details = decode_token(access_token)

        response = jsonify({
            "msg": "refresh successful",
            "access_token": access_token,
            "access_token_lifetime": access_token_details['exp'] - int(time.time()),
        })

        set_access_cookies(response, access_token)
        return response


def digby_protected():
    def wrapper(fn):
        @wraps(fn)
        def decorator(*args, **kwargs):
            try:
                verify_jwt_in_request(optional=True)
            except:
                return "Unauthorized", 403
            current_identity = get_jwt_identity()
            if current_identity or (app.config['JWT_USER'] == '' and app.config['JWT_PASSWORD'] == ''):
                return fn(*args, **kwargs)
            else:
                return "Unauthorized", 403
        return decorator
    return wrapper


