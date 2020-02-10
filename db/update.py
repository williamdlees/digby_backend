import db.salmon

def db_update():
    status = ""
    status += db.salmon.update()

    return status
