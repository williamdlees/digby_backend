from db.genomic_db import RefSeq, Feature, Sequence, SequenceFeature, Details
from sqlalchemy import or_
import os

from db.genomic_db_functions import add_feature_to_ref


# Build a GFF file for the reference sequence - contains gene-level coding and non-coding features

def build_ref_seq_gff(session, dataset_dir, ref_seq, name_prefix):
    # Feature file showing genes only: GFF3
    with open(os.path.join(dataset_dir, 'samples', name_prefix + '.gff3'), 'w') as fo:
        fo.write('##gff-version 3\n')

        #           ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
        #           ctg123 . TF_binding_site 1000  1012  .  +  .  ID=tfbs00001;Parent=gene00001

        features = session.query(Feature).filter(Feature.refseq == ref_seq)\
            .filter(Feature.feature_level == 'gene')\
            .order_by(Feature.start)\
            .all()

        for feature in features:
            if feature.feature == 'CDS' or feature.feature == 'gene':
                fo.write('%s\t.\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (ref_seq.name, 'mRNA', feature.start, feature.end, feature.strand, feature.attribute))


def build_gff(session, dataset_dir):
    details = session.query(Details).one_or_none()
    species = details.species

    ref_seqs = session.query(RefSeq).all()
    for ref_seq in ref_seqs:
        name_prefix = f"{species.replace(' ', '_')}_{ref_seq.name}"
        build_ref_seq_gff(session, dataset_dir, ref_seq, name_prefix)
        # Alignment file of all alleles and other annotated regions within samples aligned to this reference (SAM, needs external conversion to BAM)
        # FIX - add non coding regions
        with open(os.path.join(dataset_dir, 'samples', name_prefix + '.sam'), 'w') as fo:
            fo.write('@HD\tVN:1.3\tSO:coordinate\n')
            fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

            features = session.query(Feature).filter(Feature.refseq == ref_seq).order_by(Feature.start, Feature.name).all()
            for feature in features:
                if feature.feature_level == 'allele':
                    for sequence in feature.sequences:
                        legend = ('*' + sequence.name.split('*')[1] if '*' in sequence.name else sequence.name)
                        if sequence.functional == 'Functional' or not sequence.functional:
                            legend += f" ({sequence.appearances})"
                        else:
                            legend += f" ({sequence.appearances}, {sequence.functional[0]})"

                        if '.' not in sequence.sequence:
                            cigar_string = f"{len(sequence.sequence)}M"
                        else:
                            counting_d = False
                            count = 0
                            cigar_string = ''
                            for i in range(len(sequence.sequence)):
                                if counting_d:
                                    if sequence.sequence[i] == '.':
                                        count += 1
                                    else:
                                        if count:
                                            cigar_string += f"{count}D"
                                        counting_d = False
                                        count = 1
                                else:
                                    if sequence.sequence[i] != '.':
                                        count += 1
                                    else:
                                        if count:
                                            cigar_string += f"{count}M"
                                        counting_d = True
                                        count = 1

                            if count > 0:
                                if counting_d:
                                    cigar_string += f"{count}D"
                                else:
                                    cigar_string += f"{count}M"

                            sequence.sequence = sequence.sequence.replace('.', '')
                            session.commit()

                        if len(sequence.sequence) > 0:
                            fo.write('%s\t0\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t*\tNM:Z:%s\n' % (sequence.name, ref_seq.name, feature.start, cigar_string, sequence.sequence, legend))

        # Alignment file of all alleles within IMGT, plus novel alleles from samples aligned to this reference (SAM, needs external conversion to BAM)
        with open(os.path.join(dataset_dir, 'samples', name_prefix + '_imgt.sam'), 'w') as fo:
            fo.write('@HD\tVN:1.3\tSO:coordinate\n')
            fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

            features = session.query(Feature).filter(Feature.refseq == ref_seq).filter(Feature.attribute.like('%REGION%')).order_by(Feature.start).all()

            order = 1
            for feature in features:
                for sequence in feature.sequences:
                    if sequence.novel:
                        if sequence.type in ('V-REGION', 'D-REGION', 'J-REGION'):
                            gene_name = sequence.name.split('*')[1] if '*' in sequence.name else sequence.name
                            if len(sequence.sequence) > 0:
                                fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\tNM:Z:%s\n' %
                                         (sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence, order, 'novel *' + gene_name))
                        else:
                            if len(sequence.sequence) > 0:
                                fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' %
                                     (sequence.name, ref_seq.name, feature.start, len(sequence.sequence), sequence.sequence, order))
                        order += 1

                # TODO need a better way of identifying alleles here. Can't really rely on syntax. I think we need to store
                # the root name and allele as separate fields in the sequence object so that we can cope with different syntax
                imgts = session.query(Sequence)\
                    .join(SequenceFeature, SequenceFeature.sequence_id == Sequence.id) \
                    .join(Feature, Feature.id == SequenceFeature.feature_id) \
                    .join(RefSeq).filter(RefSeq.name == ref_seq.name)\
                    .filter(Feature.refseq_id == RefSeq.id) \
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
                        fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\tNM:Z:%s\n' % (imgt.name, ref_seq.name, feature.start, len(nt_sequence), nt_sequence, order, legend))
                    else:
                        fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' % (imgt.name, ref_seq.name, feature.start, len(nt_sequence), nt_sequence, order))
                    order += 1



def build_fake_ref(session, dataset_dir):
    details = session.query(Details).one_or_none()
    species = details.cell_species_label

    with open(os.path.join(dataset_dir, f'fake_{species}.fasta'), 'w') as fo:
        fo.write(f'>{species}\n')
        feats = session.query(Feature).join(RefSeq).join(SequenceFeature).join(Sequence).filter(Sequence.type.like('%REGION')).filter(RefSeq.name == 'Human_IGH').order_by(Feature.start).all()

        i = 1
        for feat in feats:
            if i <= feat.start and feat.sequences:
                seq_01 = session.query(Sequence).filter(Sequence.name == feat.name.split('_')[0] + '*01').one_or_none()
                if seq_01:
                    fo.write('n'*(feat.start-i))
                    i += (feat.start-i)
                    fo.write(seq_01.sequence)
                    i += len(seq_01.sequence)
                else:
                    seqs = session.query(Sequence).filter(Sequence.name.like(feat.name + '*%')).all()
                    if seqs:
                        fo.write('n'*(feat.start-i))
                        i += (feat.start-i)
                        fo.write(seqs[0].sequence)
                        i += len(seqs[0].sequence)
                    else:
                        print('Human reference allele %s not found.' % (feat.name + '*01'))

