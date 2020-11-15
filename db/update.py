import db.salmon
import db.igenotyper
import db.imgt_ref

def db_update():
    db.imgt_ref.update_imgt()

    status = ""
    status += db.salmon.update()
    # status += db.igenotyper.update()

    return status



