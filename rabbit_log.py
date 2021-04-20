# Copyright William Lees
#
# Very simple handler to reflect log messages onto rabbitMQ


import logging
import pika


class FlaskRabbitLogHandler(logging.Handler):
    def __init__(self, url, appname, *args, **kwargs):
        super(FlaskRabbitLogHandler, self).__init__(*args, **kwargs)
        try:
            self.url = url
            self.appname = appname
            self.connection = pika.BlockingConnection(
                    pika.URLParameters(self.url))
            self.channel = self.connection.channel()
            self.channel.queue_declare(queue='mail')
        except:
            print('Exception initialising RabbitMQ')

    def __del__(self):
        if self.connection:
            self.connection.close()

    def emit(self, record):
        try:
            self.channel.basic_publish(
                exchange='',
                routing_key='mail',
                body=self.format(record),
                properties=pika.BasicProperties(type=record.levelname, app_id=self.appname)
            )
        except:
            print('Exception logging error to RabbitMQ')
