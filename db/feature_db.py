from app import sql_db as db


class Species(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))


class RefSeq(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    locus = db.Column(db.String(100))
    sequence = db.Column(db.Text(10000000))
    length = db.Column(db.Integer)
    chromosome = db.Column(db.String(10))
    start = db.Column(db.BigInteger)
    end = db.Column(db.BigInteger)
    reference = db.Column(db.String(500))
    species_id = db.Column(db.Integer, db.ForeignKey('species.id'))
    species = db.relationship('Species', backref='refseqs')


class DataSet(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    locus = db.Column(db.String(100))
    species_id = db.Column(db.Integer, db.ForeignKey('species.id'))
    species = db.relationship('Species', backref='datasets')


class SampleSequence(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    sample_id = db.Column(db.Integer, db.ForeignKey('sample.id'), nullable=False)
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequence.id'), nullable=False)
    sample = db.relationship('Sample', backref='sequence_associations', cascade="all")
    sequence = db.relationship('Sequence', backref='sample_associations', cascade="all")
    chromosome = db.Column(db.String(50))
    chromo_count = db.Column(db.Integer)


class SequenceFeature(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    feature_id = db.Column(db.Integer, db.ForeignKey('feature.id'), nullable=False)
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequence.id'), nullable=False)
    feature = db.relationship('Feature', backref='sequence_associations', cascade="all")
    sequence = db.relationship('Sequence', backref='feature_associations', cascade="all")


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
    sequences = db.relationship('Sequence', secondary='sequence_feature')


class Sequence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    imgt_name = db.Column(db.String(100))
    type = db.Column(db.String(100))
    novel = db.Column(db.Boolean)
    deleted = db.Column(db.Boolean)
    functional = db.Column(db.String(1))
    sequence = db.Column(db.Text(10000000))
    gapped_sequence = db.Column(db.Text(10000000))
    species_id = db.Column(db.Integer, db.ForeignKey('species.id'))
    species = db.relationship('Species', backref='sequences')
    samples = db.relationship('Sample', secondary='sample_sequence')
    features = db.relationship('Feature', secondary='sequence_feature')


class Sample(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    type = db.Column(db.String(100))
    description = db.Column(db.String(500))
    date = db.Column(db.DateTime, nullable=False)
    study_id = db.Column(db.Integer, db.ForeignKey('study.id'))
    sequences = db.relationship("Sequence", secondary='sample_sequence')
    species_id = db.Column(db.Integer, db.ForeignKey('species.id'))
    ref_seq_id = db.Column(db.Integer, db.ForeignKey('ref_seq.id'))
    data_set_id = db.Column(db.Integer, db.ForeignKey('data_set.id'))
    report_link = db.Column(db.String(200))
    annot_method = db.Column(db.String(200))
    annot_ref = db.Column(db.String(500))
    data_set = db.relationship("DataSet", backref='samples')


class Study(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False)
    description = db.Column(db.String(500))
    institute = db.Column(db.String(500))
    researcher = db.Column(db.String(200))
    reference = db.Column(db.String(500))
    contact = db.Column(db.String(200))
    samples = db.relationship("Sample", backref='study')




