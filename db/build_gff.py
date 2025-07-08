from db.genomic_db import RefSeq, Feature, Sequence, SequenceFeature, Details
from sqlalchemy import or_, and_
import os
from Bio import Align

from receptor_utils import simple_bio_seq as simple
from db.cigar import Cigar


def aligned_diff(novel_seq: str, ref_seq: str):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.match_score = 2
    aligner.mismatch_score = 0
    aligner.target_end_gap_score = 1

    best_alignment = aligner.align(novel_seq, ref_seq)[0]

    # return 0-based start, end and score
    return best_alignment[0], best_alignment[1], best_alignment.score


def create_cigar_from_alignment(alignment):
    """
    Create a CIGAR string from a pairwise alignment.

    Args:
        alignment: A tuple with (query_aligned, reference_aligned) from aligned_diff function

    Returns:
        A string representing the CIGAR operations
    """
    query_aligned, ref_aligned = alignment[0], alignment[1]

    cigar_ops = []
    current_op = None
    count = 0

    # Walk through the alignment paths
    for q_char, r_char in zip(query_aligned, ref_aligned):
        # Match or mismatch (both are 'M' in CIGAR)
        if q_char != '-' and r_char != '-':
            op = 'M'
        # Insertion in query (gap in reference)
        elif q_char != '-' and r_char == '-':
            op = 'I'
        # Deletion in query (gap in query)
        elif q_char == '-' and r_char != '-':
            op = 'D'
        else:
            # This shouldn't happen in normal alignments
            continue

        # If we're continuing the same operation, increment the count
        if op == current_op:
            count += 1
        else:
            # If we're changing operations, save the previous one
            if current_op:
                cigar_ops.append(f"{count}{current_op}")
            # Start counting the new operation
            current_op = op
            count = 1

    # Add the final operation
    if current_op:
        cigar_ops.append(f"{count}{current_op}")

    # Join all operations to create the CIGAR string
    cigar_string = "".join(cigar_ops)

    return cigar_string


def build_gff(session, dataset_dir):
    details = session.query(Details).one_or_none()
    species = details.species

    ref_seqs = session.query(RefSeq).all()
    for ref_seq in ref_seqs:
        name_prefix = f"{species.replace(' ', '_')}_{ref_seq.name.split(':')[0]}"
        build_ref_seq_gff(session, dataset_dir, ref_seq, name_prefix)
        phased_feature_alignment(dataset_dir, name_prefix + '_phased', ref_seq, session)
        unphased_feature_alignment(dataset_dir, name_prefix, ref_seq, session)
        all_imgt_and_novel_v_region_alignment(dataset_dir, name_prefix, ref_seq, session)


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


# Alignment file of all alleles within IMGT, plus novel alleles from samples aligned to this reference (SAM, needs external conversion to BAM)
# TODO - this hasn't been used for a while and is probably broken
def all_imgt_and_novel_v_region_alignment(dataset_dir, name_prefix, ref_seq, session):
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
            imgts = session.query(Sequence) \
                .join(SequenceFeature, SequenceFeature.sequence_id == Sequence.id) \
                .join(Feature, Feature.id == SequenceFeature.feature_id) \
                .join(RefSeq).filter(RefSeq.name == ref_seq.name) \
                .filter(Feature.refseq_id == RefSeq.id) \
                .filter(or_(Sequence.name.like(feature.name.split('_')[0] + '*%'), Sequence.name == feature.name)) \
                .filter(Sequence.novel == False).order_by(Sequence.name.desc()).all()

            for imgt in imgts:
                nt_sequence = imgt.sequence

                if imgt.gapped_sequence and imgt.gapped_sequence[0] == '.' and imgt.type == 'V-REGION':
                    i = 0
                    pad = ''
                    while imgt.gapped_sequence[i] == '.':
                        if imgts[-1].gapped_sequence[i] != '.':  # let's just assume that *01 is never truncated
                            pad += 'x'
                        i += 1
                    nt_sequence = pad + nt_sequence

                if imgt.type in ('V-REGION', 'D-REGION', 'J-REGION'):
                    legend = ('*' + imgt.name.split('*')[1]) if '*' in imgt.name else imgt.name
                    fo.write('%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\tNM:Z:%s\n' % (
                        imgt.name, ref_seq.name, feature.start, len(nt_sequence), nt_sequence, order, legend))
                else:
                    fo.write(
                        '%s\t0\t%s\t%d\t255\t%dM\t*\t0\t0\t%s\t*\tOD:i:%d\n' % (imgt.name, ref_seq.name, feature.start, len(nt_sequence), nt_sequence, order))
                order += 1


# Alignment file of all alleles and other annotated regions within samples aligned to this reference (SAM, needs external conversion to BAM)
# Unphased - just push out all sequences we see for each feature
def unphased_feature_alignment(dataset_dir, name_prefix, ref_seq, session):
    with open(os.path.join(dataset_dir, 'samples', name_prefix + '.sam'), 'w') as fo:
        fo.write('@HD\tVN:1.3\tSO:coordinate\n')
        fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

        features = session.query(Feature).filter(Feature.refseq == ref_seq).order_by(Feature.start, Feature.name).all()
        for feature in features:
            if feature.feature_level == 'allele' and feature.feature_type != 'gene_sequence':
                for sequence in feature.sequences:
                    if sequence.type == 'C-REGION':
                        continue
                    feature_name = ('*' + sequence.name.split('*')[1].replace('_phased', '') if '*' in sequence.name else sequence.name)
                    rec = feature_gff_rec(feature, feature_name, ref_seq, sequence, session)
                    if rec:
                        fo.write(rec)


# Alignment file of all alleles and other annotated regions within samples aligned to this reference (SAM, needs external conversion to BAM)
# Phased - report each combination of features/sequences we see, with subject counts
def phased_feature_alignment(dataset_dir, name_prefix, ref_seq, session):
    with open(os.path.join(dataset_dir, 'samples', name_prefix + '.sam'), 'w') as fo:
        fo.write('@HD\tVN:1.3\tSO:coordinate\n')
        fo.write('@SQ\tSN:%s\tLN:%d\n' % (ref_seq.name, len(ref_seq.sequence)))

        features = session.query(Feature).filter(and_(Feature.refseq == ref_seq, Feature.feature_type == 'gene_sequence')).order_by(Feature.start, Feature.name).all()

        if len(features) == 0:
            print('gene_sequence features not found. falling back to REGIONS')
            features = session.query(Feature).filter(and_(Feature.refseq == ref_seq, Feature.feature_type.like('%REGION'))).order_by(Feature.start, Feature.name).all()

        for feature in features:
            if feature.feature_level == 'allele':
                for sequence in feature.sequences:
                    feature_name = '*' + sequence.name.split('*')[1]
                    rec = feature_gff_rec(feature, feature_name, ref_seq, sequence, session)
                    if rec:
                        fo.write(rec)


def feature_gff_rec(feature, feature_name, ref_seq, sequence, session):
    legend = feature_name
    if sequence.functional == 'Functional' or not sequence.functional:
        legend += f" ({sequence.appearances})"
    else:
        legend += f" ({sequence.appearances}, {sequence.functional[0]})"
    # IGenotyper records will always have a cigar string
    
    cigar_string = feature.feature_cigar
    
    if len(sequence.sequence) > 0 and cigar_string and cigar_string != '':
        if 'D' in cigar_string:
            sequence.sequence = sequence.sequence.replace('-', '')
        
        seq = sequence.sequence
        cigar_string = Cigar(cigar_string)

        if len(seq) != cigar_string.__len__():
            print(f'Warning: sequence length {len(seq)} does not match cigar length {cigar_string.__len__()} for feature {feature.name} in refseq {ref_seq.name}. Probable length mismatch between reference sequence and bed file coords')
            alignment = aligned_diff(seq, feature.feature_seq)
            cigar_string = create_cigar_from_alignment(alignment)
            
        if feature.strand != ref_seq.sense:
            seq = simple.reverse_complement(seq)

        if ref_seq.sense == '-':
            cigar_string = cigar_string._reverse_cigar()

        return '%s\t0\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t*\tNM:Z:%s\n' % (
            sequence.name, ref_seq.name, feature.start, cigar_string, seq, legend)
    else:
        return None


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

