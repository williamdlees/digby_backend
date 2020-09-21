# Flask / Celery integration
# from https://stackoverflow.com/questions/12044776/how-to-use-flask-sqlalchemy-in-a-celery-task

import flask
from flask_sqlalchemy import SQLAlchemy
from celery import Celery

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


celery = FlaskCelery('tasks', broker='pyamqp://guest@localhost//', backend='rpc://')
db = SQLAlchemy()

from db.feature_db import Species, RefSeq, Feature, Sequence, SequenceFeature, Sample, Study, SampleSequence

@celery.task()
def add_together(a, b):
    sp = db.session.query(Species).all()

    if sp:
        return [row.name for row in sp]
    else:
        return 'None'
