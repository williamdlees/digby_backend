from app import db
from db.genomic_db import RefSeq, Species, Sample, Sequence, Feature, SampleSequence, SequenceFeature, Study, DataSet


def delete_dependencies(species):
    db.session.commit()

    if species:
        sequence_features = db.session.query(SequenceFeature.id).join(Sequence).join(Species).filter(Species.name == species).all()
        delete_q = SequenceFeature.__table__.delete().where(SequenceFeature.id.in_(sequence_features))
        db.session.execute(delete_q)
        db.session.commit()

        sample_sequences = db.session.query(SampleSequence.id).join(Sequence).join(Species).filter(Species.name == species).all()
        delete_q = SampleSequence.__table__.delete().where(SampleSequence.id.in_(sample_sequences))
        db.session.execute(delete_q)
        db.session.commit()

        sequences = db.session.query(Sequence.id).join(Species).filter(Species.name == species).all()
        delete_q = Sequence.__table__.delete().where(Sequence.id.in_(sequences))
        db.session.execute(delete_q)
        db.session.commit()

        features = db.session.query(Feature.id).join(RefSeq).join(Species).filter(Species.name == species).all()
        delete_q = Feature.__table__.delete().where(Feature.id.in_(features))
        db.session.execute(delete_q)
        db.session.commit()

        studies = db.session.query(Study.id).join(Sample).join(Species).filter(Species.name == species).all()

        samples = db.session.query(Sample.id).join(Species).filter(Species.name == species).all()
        delete_q = Sample.__table__.delete().where(Sample.id.in_(samples))
        db.session.execute(delete_q)
        db.session.commit()

        delete_q = Study.__table__.delete().where(Study.id.in_(studies))
        db.session.execute(delete_q)
        db.session.commit()

        ref_seqs = db.session.query(RefSeq.id).join(Species).filter(Species.name == species).all()
        delete_q = RefSeq.__table__.delete().where(RefSeq.id.in_(ref_seqs))
        db.session.execute(delete_q)
        db.session.commit()

        data_sets = db.session.query(DataSet.id).join(Species).filter(Species.name == species).all()
        delete_q = DataSet.__table__.delete().where(DataSet.id.in_(data_sets))
        db.session.execute(delete_q)
        db.session.commit()

    else:
        SequenceFeature.query.delete()
        SampleSequence.query.delete()
        Sequence.query.delete()
        Feature.query.delete()
        Sample.query.delete()
        RefSeq.query.delete()
        Study.query.delete()

    db.session.commit()
