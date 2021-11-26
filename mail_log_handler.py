# Copyright William Lees
#
# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2, 
# the English version of which is available here: https://perma.cc/DK5U-NDVE
#

import logging

from flask_mail import Message
from flask.globals import current_app, _app_ctx_stack
from flask import Flask

class FlaskMailLogHandler(logging.Handler):

    def __init__(self, mail, sender, recipients, subject, *args, **kwargs):
        super(FlaskMailLogHandler, self).__init__(*args, **kwargs)
        self.mail = mail
        self.sender = sender
        self.recipients = recipients
        self.subject = subject

    def emit(self, record):
        if _app_ctx_stack.top is None:
            # we are outside the application context. Need to build one to send the mail
            app = Flask(__name__)
        else:
            app = current_app
        with app.app_context():
            self.mail.send(
                Message(
                    sender=self.sender,
                    recipients=self.recipients,
                    body=self.format(record),
                    subject=self.subject
                )
            )