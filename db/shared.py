from app import db
from db.feature_db import RefSeq, Species

def delete_dependencies(species):
    if species:
        ref_seqs = db.session.query(RefSeq).join(Species).filter(Species.name == species).all()
    else:
        ref_seqs = db.session.query(RefSeq).join(Species).all()

    for ref_seq in ref_seqs:
        for feature in ref_seq.features:
            if feature.sequence and feature.sequence.novel == 1 and len(feature.sequence.features) == 1:
                db.session.delete(feature.sequence)
                db.session.commit()
            feature.samples = []
            db.session.commit()
            db.session.delete(feature)
        db.session.delete(ref_seq)
    db.session.commit()
