# Mixin class defining DB property fields

from sqlalchemy import Column, Integer, String, DateTime
import datetime


class DB_PropertyMixin(object):
    id = Column(Integer, primary_key=True)
    created_on = Column(DateTime, default=datetime.datetime.now)
    created_by = Column(String(100))
    software_version = Column(String(100))
