# Copyright William Lees
#
# Very simple handler to reflect log messages onto rabbitMQ


import logging
import pika


class FlaskRabbitLogHandler(logging.Handler):
    def __init__(self, url, appname, *args, **kwargs):
        super(FlaskRabbitLogHandler, self).__init__(*args, **kwargs)
        self.url = url
        self.appname = appname

    def __del__(self):
        pass

    def emit(self, record):
        # set up the connection here, because we only expect the occasional message and
        # blocking connections are bad at maintaining the heartbeat
        self.connection = pika.BlockingConnection(
                pika.URLParameters(self.url))
        self.channel = self.connection.channel()
        self.channel.queue_declare(queue='mail')

        self.channel.basic_publish(
            exchange='',
            routing_key='mail',
            body=self.format(record),
            properties=pika.BasicProperties(type=record.levelname, app_id=self.appname)
        )

        if self.connection:
            self.connection.close()
