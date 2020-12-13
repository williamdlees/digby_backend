# Flask / Celery integration
# from https://stackoverflow.com/questions/12044776/how-to-use-flask-sqlalchemy-in-a-celery-task
import importlib
import traceback

import flask
from flask_sqlalchemy import SQLAlchemy
from celery import Celery, current_task
from werkzeug.exceptions import BadRequest


class FlaskCelery(Celery):

    def __init__(self, *args, **kwargs):

        super(FlaskCelery, self).__init__(*args, **kwargs)
        self.patch_task()

        if 'app' in kwargs:
            self.init_app(kwargs['app'])

    def patch_task(self):
        TaskBase = self.Task
        _celery = self

        class ContextTask(TaskBase):
            abstract = True

            def __call__(self, *args, **kwargs):
                if flask.has_app_context():
                    return TaskBase.__call__(self, *args, **kwargs)
                else:
                    with _celery.app.app_context():
                        return TaskBase.__call__(self, *args, **kwargs)

        self.Task = ContextTask

    def init_app(self, app):
        self.app = app
        self.config_from_object(app.config)
        return SQLAlchemy(app)



celery = FlaskCelery('tasks', broker='pyamqp://guest@localhost//', backend='redis://localhost:6379/0')


@celery.task(bind=True)
def run_report(self, report_name, format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    runner = importlib.import_module('api.reports.' + report_name)

    try:
        self.update_state(state='PENDING', meta={'stage': 'preparing data'})
        return runner.run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params)
    except BadRequest as bad:
        print('BadRequest raised during report processing: %s' % bad.description)
        return {'status': 'error', 'description': bad.description}
    except Exception as e:
        print('Exception raised during report processing: %s' % traceback.format_exc())
        return {'status': 'error', 'description': 'Unexpected error when running report: %s' % traceback.format_exc()}

