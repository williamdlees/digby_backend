from app import db
from db.feature_db import RefSeq, Species, Sample, Sequence, Feature, SampleSequence, SequenceFeature, Study

def delete_dependencies(species):
    db.session.commit()

    if species:
        ref_seqs = db.session.query(RefSeq).join(Species).filter(Species.name == species).all()

        for ref_seq in ref_seqs:
            for feature in ref_seq.features:
                feature.sequences = []
                feature.samples = []
            db.session.commit()
            for feature in ref_seq.features:
                db.session.delete(feature)
            db.session.commit()
            db.session.delete(ref_seq)

    else:
        SampleSequence.query.delete()
        SequenceFeature.query.delete()
        Sequence.query.delete()
        Feature.query.delete()
        Sample.query.delete()
        RefSeq.query.delete()
        Study.query.delete()


    db.session.commit()
