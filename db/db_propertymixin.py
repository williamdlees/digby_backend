# Mixin class defining DB property fields

from sqlalchemy import Column, Integer, String, DateTime
import datetime


class Details_Mixin(object):
    id = Column(Integer, primary_key=True)
    dbtype = Column(String(100))
    species = Column(String(100))
    locus = Column(String(100))
    created_on = Column(DateTime, default=datetime.datetime.now)
    created_by = Column(String(100))
    software_commit_id = Column(String(100))
    software_branch = Column(String(100))
