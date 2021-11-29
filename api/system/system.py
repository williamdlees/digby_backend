# Services related to vdjbase repseq-based data sets
import requests
from flask import jsonify
from flask_jwt_extended import create_access_token, set_access_cookies
from flask_restx import Resource, reqparse, fields, marshal, inputs
from api.restx import api
import json
from app import vdjbase_dbs, app, db


ns = api.namespace('system', description='System-specific information')


@ns.route('/config')
class ConfigApi(Resource):
    def get(self):
        """ Return internal configuration details """

        config = {'wp_news_url': app.config['WORDPRESS_NEWS_URL'], 'wp_help_url': app.config['WORDPRESS_HELP_URL'], 'wp_rest': app.config['WORDPRESS_REST']}

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
            access_token = create_access_token(identity="example_user")
            response = jsonify({"msg": "login successful", "token": access_token})
            set_access_cookies(response, access_token)
        else:
            response = jsonify({"msg": "invalid username or password"}), 401
        return response
