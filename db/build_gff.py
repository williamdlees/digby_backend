from app import db
from db.feature_db import RefSeq, Feature, Species, Sequence, SequenceFeature


def build_gffs():
    build_gff('Atlantic Salmon')
    build_gff('Human')
    build_fake_human_ref()
    return('GFFs built!')

def build_gff(species):
    ref_seqs = db.session.query(RefSeq).join(Species).filter(Species.name == species).all()

    for ref_seq in ref_seqs:
        # Feature file showing genes only: GFF3
        with open('static/gff/%s_%s.gff3' % (species.replace(' ', '_'), ref_seq.name), 'w') as fo:
            fo.write('##gff-version 3\n')

#           ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
#           ctg123 . TF_binding_site 1000  1012  .  +  .  ID=tfbs00001;Parent=gene00001

            features = db.session.query(Feature).filter(Feature.refseq == ref_seq).order_by(Feature.start).all()
            for feature in features:
                if feature.feature == 'gene':
                    fo.write('%s\t.\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (ref_seq.name, 'mRNA', feature.start, feature.end, feature.strand, feature.attribute))

        # Alignment file of all alleles and other annotated regions within samples aligned to this reference (SAM, needs external conversion to BAM)
        # FIX - add non coding regions
        with open('static/gff/%s_%s.sam' % (species.replace(' ', '_'), ref_seq.name), 'w') as fo:
            fo.write('@HD\tVN:1.3\tSO:coordinate\n')
            fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

            for feature in features:
                if feature.feature != 'gene':
                    for sequence in feature.sequences:
                        fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\n' %(sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence))

        # Alignment file of all alleles within IMGT, plus novel alleles from samples aligned to this reference (SAM, needs external conversion to BAM)
        with open('static/gff/%s_%s_imgt.sam' % (species.replace(' ', '_'), ref_seq.name), 'w') as fo:
            fo.write('@HD\tVN:1.3\tSO:coordinate\n')
            fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

            features = db.session.query(Feature).filter(Feature.refseq == ref_seq).filter(Feature.name.like('%REGION')).order_by(Feature.start).all()

            for feature in features:
                imgts = db.session.query(Sequence).join(Species).join(RefSeq).filter(RefSeq.name == ref_seq.name).\
                    filter(Sequence.name.like(feature.name.split('_')[0] + '*%')).filter(Sequence.novel == False).order_by(Sequence.name).all()

                order = 1
                for imgt in imgts:
                    fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' %(imgt.name, ref_seq.name, feature.start, len(imgt.sequence), imgt.sequence, order))
                    order += 1

                order = 100
                for sequence in feature.sequences:
                    if sequence.novel:
                        fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' %(sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence, order))


def build_fake_human_ref():
    with open('static/gff/Human_IMGT_IGH.fasta', 'w') as fo:
        fo.write('>IMGT_IGH\n')
        feats = db.session.query(Feature).join(RefSeq).join(SequenceFeature).join(Sequence).filter(Sequence.type.like('%REGION')).filter(RefSeq.name == 'Human_IGH').order_by(Feature.start).all()

        i = 1
        for feat in feats:
            if i <= feat.start and feat.sequences:
                seq_01 = db.session.query(Sequence).join(Species).filter(Species.name == 'Human').filter(Sequence.name == feat.name.split('_')[0] + '*01').one_or_none()
                if seq_01:
                    fo.write('n'*(feat.start-i))
                    i += (feat.start-i)
                    fo.write(seq_01.sequence)
                    i += len(seq_01.sequence)
                else:
                    seqs = db.session.query(Sequence).filter(Sequence.name.like(feat.name + '*%')).all()
                    if seqs:
                        fo.write('n'*(feat.start-i))
                        i += (feat.start-i)
                        fo.write(seqs[0].sequence)
                        i += len(seqs[0].sequence)
                    else:
                        print('Human reference allele %s not found.' % (feat.name + '*01'))

