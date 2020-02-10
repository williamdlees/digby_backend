

from app import db

class Species(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))


class RefSeq(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    locus = db.Column(db.String(100))
    sequence = db.Column(db.Text(10000000))
    length = db.Column(db.Integer)
    species_id = db.Column(db.Integer, db.ForeignKey('species.id'))
    species = db.relationship('Species', backref='refseqs')


class Feature(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    feature = db.Column(db.String(100))
    start = db.Column(db.Integer)
    end = db.Column(db.Integer)
    score = db.Column(db.Float)
    strand = db.Column(db.String(1))
    frame = db.Column(db.Integer)
    attribute = db.Column(db.String(1000))
    feature_id = db.Column(db.Integer)
    parent_id = db.Column(db.Integer)
    refseq_id = db.Column(db.Integer, db.ForeignKey('ref_seq.id'))
    refseq = db.relationship('RefSeq', backref='features')
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequence.id'))
    sequence = db.relationship('Sequence', backref='features')


class Sequence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    imgt_name = db.Column(db.String(100))
    type = db.Column(db.String(100))
    sequence = db.Column(db.Text(10000000))
    gapped_sequence = db.Column(db.Text(10000000))
    refseq_id = db.Column(db.Integer, db.ForeignKey('ref_seq.id'))
    refseq = db.relationship('RefSeq', backref='sequences')
