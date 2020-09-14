# Services related to vdjbase repseq-based data sets

from flask import request, redirect, send_file
from flask_restplus import Resource, reqparse, fields, marshal, inputs
from api.restplus import api
import json
from werkzeug.exceptions import BadRequest
from db.feature_db import RefSeq, Sample
from db.vdjbase_model import Sample as vdjb_Sample
from api.genomic.genomic import find_genomic_samples, find_genomic_filter_params
from api.vdjbase.vdjbase import find_vdjbase_samples, find_rep_filter_params
from app import app
import importlib
import subprocess
import os
import traceback
from json.decoder import JSONDecodeError
import tempfile

SYSDATA = os.path.join(app.config['R_SCRIPT_PATH'], 'sysdata.rda')

ns = api.namespace('reports', description='Reports that can be run on subsets of genomic and repseq data')
report_defs = None

def load_report_defs():
    with open('api/reports/reports.json') as fi:
        global report_defs
        report_defs = json.load(fi)['reports']
        for k, report_def in report_defs.items():
            if 'thumbnail' in report_def:
                report_defs[k]['thumbnail'] = app.config['STATIC_LINK'] + 'img/reports/' + report_def['thumbnail']

report_list_arguments = reqparse.RequestParser()
report_list_arguments.add_argument('species', type=str)
report_list_arguments.add_argument('genomic_datasets', type=str)
report_list_arguments.add_argument('rep_datasets', type=str)

@ns.route('/reports/list')
@api.response(404, 'No reports available!')
class ReportsApi(Resource):
    @api.expect(report_list_arguments, validate=True)
    def get(self):
        try:
            args = report_list_arguments.parse_args(request)
            scope = set()

            if len(args.genomic_datasets):
                genomic_datasets = args.genomic_datasets.split(',')
                scope |= {'gen_sample'}
            else:
                genomic_datasets = None

            if len(args.rep_datasets):
                rep_datasets = args.rep_datasets.split(',')
                scope |= {'rep_sample'}
            else:
                rep_datasets = None

            genomic_filter_params = find_genomic_filter_params(args.species, genomic_datasets) if genomic_datasets is not None else []
            rep_filter_params = find_rep_filter_params(args.species, rep_datasets) if rep_datasets is not None else []

            available_reports = {}

            for k,v in report_defs.items():
                if len(set(v['scope']) & scope):
                    available_reports[k] = v

            combined_filter_params = {}

            for params_list in (genomic_filter_params, rep_filter_params):
                if params_list is not None:
                    for param in params_list:
                        if param['id'] not in combined_filter_params:
                            combined_filter_params[param['id']] = param
                        elif 'options' in param:
                            combined_filter_params[param['id']]['options'] = list(set(combined_filter_params[param['id']]['options'] + param['options']))

            return {
                'reports': available_reports,
                'filters': {
                    'combined': combined_filter_params,
                    'gen': genomic_filter_params,
                    'rep': rep_filter_params
                }
            }

        except BadRequest as bad:
            print('BadRequest raised during report processing: %s' % bad.description)
            raise bad
        except Exception:
            print('Exception encountered processing report request: %s' % traceback.format_exc())
            raise BadRequest('Error encountered while processing report request')


report_arguments = reqparse.RequestParser()
report_arguments.add_argument('format', type=str)           # pdf or html
report_arguments.add_argument('species', type=str)
report_arguments.add_argument('genomic_datasets', type=str)
report_arguments.add_argument('genomic_filters', type=str)
report_arguments.add_argument('rep_datasets', type=str)
report_arguments.add_argument('rep_filters', type=str)
report_arguments.add_argument('params', type=str)

@ns.route('/reports/run/<string:report_name>')
@api.response(404, 'Malformed request')
class ReportsRunApi(Resource):
    @api.expect(report_arguments, validate=True)
    def get(self, report_name):
        try:
            args = report_arguments.parse_args(request)

            if report_name not in report_defs:
                raise BadRequest('No such report')

            try:
                genomic_datasets = args.genomic_datasets.split(',') if len(args.genomic_datasets) else None
                genomic_filters = json.loads(args.genomic_filters)

                if genomic_datasets is not None:
                    genomic_samples = find_genomic_samples([Sample.name, RefSeq.name], args.species, genomic_datasets, genomic_filters)
                else:
                    genomic_samples = []

                rep_datasets = args.rep_datasets.split(',') if len(args.rep_datasets) else None
                rep_filters = json.loads(args.rep_filters)

                if rep_datasets is not None:
                    rep_samples = find_vdjbase_samples([vdjb_Sample.name, vdjb_Sample.id], args.species, rep_datasets, rep_filters)
                else:
                    rep_samples = []
            except:
                raise BadRequest("Malformed request")

            params = json.loads(args.params)

            if len(genomic_samples)== 0 and len(rep_samples) == 0:
                raise BadRequest('No samples selected')

            # maybe we should check types as well
            for p in report_defs[report_name]['params']:
                if p['id'] not in params.keys():
                    raise BadRequest('Missing parameter: %s' % p)

            runner = importlib.import_module('api.reports.' + report_name)
        except JSONDecodeError:
            print('Exception encountered processing JSON-encoded field: %s' % traceback.format_exc())
            raise BadRequest('Error encountered while processing JSON-encoded field')
        except BadRequest as bad:
            print('BadRequest raised during report processing: %s' % bad.description)
            raise bad
        except Exception:
            print('Exception encountered processing report request: %s' % traceback.format_exc())
            raise BadRequest('Error encountered while processing report request')


        try:
            return runner.run(args.format, args.species, genomic_samples, rep_samples, params)
        except BadRequest as bad:
            print('BadRequest raised during report processing: %s' % bad.description)
            raise bad
        except Exception:
            print('Exception raised during report processing: %s' % traceback.format_exc())
            raise BadRequest('Error encountered while processing report')


# Make a unique file in the output directory
def make_output_file(format):
    with tempfile.NamedTemporaryFile(suffix='.' + format, dir=app.config['OUTPUT_PATH'], delete=False) as fo:
        output_path = fo.name
        fo.close()
    return output_path


# R Script Runner
def run_rscript(script, args, cwd=app.config['R_SCRIPT_PATH']):
    cmd_line = ['Rscript', os.path.join(app.config['R_SCRIPT_PATH'], script)]
    cmd_line.extend(args)
    print("Running Rscript: '%s'\n" % ' '.join(cmd_line))
    proc = subprocess.Popen(cmd_line, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    (stdout, stderr) = proc.communicate()

    print(stdout.decode("utf-8"))

    if proc.returncode != 0:
        msg = str(stderr.decode("utf-8"))
        if "Execution halted" in msg:
            print("Got error in Rscript:\n" + msg)
            msg = msg.split('\n')
            errors = []
            for line in msg:
                if 'rror' in line:
                    errors.append(line)
            raise BadRequest('Error running report: %s' % '\n'.join(errors))

    return True


# Send a report - the file should be in OUTPUT_PATH
def send_report(filename, format, attachment_filename=None):
    return {'status': 'ok', 'filename': attachment_filename, 'url': app.config['OUTPUT_REPORT_LINK'] + os.path.basename(filename)}


