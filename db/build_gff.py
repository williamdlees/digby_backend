from app import db
from db.feature_db import RefSeq, Feature, Species, Sequence, SequenceFeature
from sqlalchemy import or_


def build_gffs():

    species = db.session.query(Species.name.distinct()).join(RefSeq).all()
    species = [x[0] for x in species]

    for sp in species:
        build_gff(sp)
    return('GFFs built! Now run /mnt/d/Research/digby_backend/static/gff/make_bam')

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
                if feature.feature == 'CDS':
                    fo.write('%s\t.\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (ref_seq.name, 'mRNA', feature.start, feature.end, feature.strand, feature.attribute))

        # Alignment file of all alleles and other annotated regions within samples aligned to this reference (SAM, needs external conversion to BAM)
        # FIX - add non coding regions
        with open('static/gff/%s_%s.sam' % (species.replace(' ', '_'), ref_seq.name), 'w') as fo:
            fo.write('@HD\tVN:1.3\tSO:coordinate\n')
            fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

            features = db.session.query(Feature).filter(Feature.refseq == ref_seq).order_by(Feature.start).all()
            for feature in features:
                if feature.feature == 'CDS':
                    for sequence in feature.sequences:
                        if sequence.type in ('V-REGION', 'D-REGION', 'J-REGION'):  # '_' is a fudge for salmon line-bred, no alleles
                            legend = 'novel ' if sequence.novel else ''
                            legend = legend + ('*' + sequence.name.split('*')[1] if '*' in sequence.name else sequence.name)
                            fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tNM:Z:%s\n' %(sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence, legend))
                        else:
                            fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\n' %(sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence))

        # Alignment file of all alleles within IMGT, plus novel alleles from samples aligned to this reference (SAM, needs external conversion to BAM)
        with open('static/gff/%s_%s_imgt.sam' % (species.replace(' ', '_'), ref_seq.name), 'w') as fo:
            fo.write('@HD\tVN:1.3\tSO:coordinate\n')
            fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

            features = db.session.query(Feature).filter(Feature.refseq == ref_seq).filter(Feature.attribute.like('%REGION%')).order_by(Feature.start).all()

            order = 1
            for feature in features:
                for sequence in feature.sequences:
                    if sequence.novel:
                        if sequence.type in ('V-REGION', 'D-REGION', 'J-REGION'):
                            gene_name = sequence.name.split('*')[1] if '*' in sequence.name else sequence.name
                            fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\tNM:Z:%s\n' %
                                     (sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence, order, 'novel *' + gene_name))
                        else:
                            fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' %
                                     (sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence, order))
                        order += 1

                # TODO need a better way of identifying alleles here. Can't really rely on syntax. I think we need to store
                # the root name and allele as separate fields in the sequence object so that we can cope with different syntax
                imgts = db.session.query(Sequence).join(Species).join(RefSeq).filter(RefSeq.name == ref_seq.name)\
                    .filter(or_(Sequence.name.like(feature.name.split('_')[0] + '*%'), Sequence.name == feature.name))\
                    .filter(Sequence.novel == False).order_by(Sequence.name.desc()).all()

                for imgt in imgts:
                    nt_sequence = imgt.sequence

                    if imgt.gapped_sequence and imgt.gapped_sequence[0] == '.' and imgt.type == 'V-REGION':
                        i = 0
                        pad = ''
                        while imgt.gapped_sequence[i] == '.':
                            if imgts[-1].gapped_sequence[i] != '.':              # let's just assume that *01 is never truncated
                                pad += 'x'
                            i += 1
                        nt_sequence = pad + nt_sequence

                    if imgt.type in ('V-REGION', 'D-REGION', 'J-REGION'):
                        legend = ('*' + imgt.name.split('*')[1]) if '*' in imgt.name else imgt.name
                        fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\tNM:Z:%s\n' %(imgt.name, ref_seq.name, feature.start, len(nt_sequence), nt_sequence, order, legend))
                    else:
                        fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' %(imgt.name, ref_seq.name, feature.start, len(nt_sequence), nt_sequence, order))
                    order += 1



def build_fake_human_ref():
    with open('static/gff/Human_Human_IGH.fasta', 'w') as fo:
        fo.write('>Human_IGH\n')
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

