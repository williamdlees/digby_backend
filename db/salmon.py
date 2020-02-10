from lxml import html
import requests

from app import db
from db.feature_db import Species, RefSeq, Feature, Sequence
from sqlalchemy.types import Text

type = {
    "5'UTR": 'five_prime_UTR',
    'L-PART1': 'CDS',
    'V-INTRON': 'intron',
    'L-PART2': 'CDS',
    'V-REGION': 'CDS',
    'D-REGION': 'CDS',
    'J-REGION': 'CDS',
    "3'UTR": 'three_prime_UTR',
}


def delete_dependencies(ctg):
    ref_seqs = db.session.query(RefSeq).filter_by(name=ctg)

    for ref_seq in ref_seqs:
        for sequence in ref_seq.sequences:
            for feature in sequence.features:
                db.session.delete(feature)
            db.session.delete(sequence)

        for feature in ref_seq.features:        # because not all features are associated with sequences
            db.session.delete(feature)

        db.session.delete(ref_seq)
        db.session.commit()


def update():
    ctg = 'GU129139'
    ret = "Importing Salmon IGHV from %s\n" % ctg

    delete_dependencies(ctg)

    page = requests.get('http://www.imgt.org/ligmdb/view?id=GU129139')
    tree = html.fromstring(page.content)

    seq_text = tree.xpath('//div[@class="sequence"]/pre')[0]

    sequence = ''
    for row in seq_text.text.split('\n'):
        if len(row) > 75:
            sequence += row[1:70].replace(' ', '')

    sp = db.session.query(Species).filter_by(name='Atlantic Salmon').one_or_none()

    if not sp:
        sp = Species(name='Atlantic Salmon')
        db.session.add(sp)

    ref_seq = RefSeq(name=ctg, locus='IGH', species=sp, sequence=sequence, length=len(sequence))
    db.session.add(ref_seq)
    db.session.commit()

    features = tree.xpath('//div[@class="features"]/table')[0]
    rows = iter(features)

    state = None
    name = None
    gene_range = None
    strand = None
    parent_id = 0

    for row in rows:
        values = [col.text for col in row]

        if len(values) < 3:
            continue

        def get_range(s):
            gene_range = s.split('..')
            if len(gene_range) < 2:
                print('Invalid gene range found: %s' % s)
                return (('1', '1'), '+')

            for i in (0, 1):
                gene_range[i] = gene_range[i].replace('>', '').replace('<', '')

            strand = '+'

            if 'complement(' in gene_range[0]:
                gene_range[0] = gene_range[0].replace('complement(', '')
                gene_range[1] = gene_range[1].replace(')', '')
                strand = '-'

            try:
                if int(gene_range[1]) - int(gene_range[0]) > 10000000 or int(gene_range[0]) > int(gene_range[1]):
                    print('Invalid gene range found: %s' % s)
                    return (('1', '1'), '+')
            except:
                print('Invalid gene range found: %s' % s)
                return (('1', '1'), '+')

            return (gene_range, strand)

        if not state and values[0] in ['V-GENE', 'D-GENE', 'J-GENE']:
            gene_range, strand = get_range(values[2])
            state = values[0]

        elif state and not name:
            if values[1] == 'IMGT_allele':
                name = values[2]
                parent_id += 1
                fp = Feature(
                    name=name,
                    feature='gene',
                    start=gene_range[0],
                    end=gene_range[1],
                    strand=strand,
                    attribute='Name=%s;ID=%s' % (name, parent_id),
                    feature_id=parent_id,
                )
                ref_seq.features.append(fp)
                parent_id += 1

                f = Feature(
                    name=name,
                    feature='mRNA',
                    start=gene_range[0],
                    end=gene_range[1],
                    strand=strand,
                    attribute='Name=%s;ID=%s;Parent=%s' % (name, parent_id, parent_id-1),
                    feature_id=parent_id,
                    parent_id=parent_id-1,
                )
                ref_seq.features.append(f)

        elif state and name:
            if (state == 'V-GENE' and values[0] in ["5'UTR", 'L-PART1', 'V-INTRON', 'L-PART2', 'V-REGION', "3'UTR"]) \
                    or (state == 'D-GENE' and values[0] in ["5'UTR", 'D-REGION', "3'UTR"]) \
                    or (state == 'J-GENE' and values[0] in ["5'UTR", 'D-REGION', "3'UTR"]):
                gene_range, strand = get_range(values[2])

                s = Sequence(
                    name=name + '_' + values[0],
                    imgt_name=name,
                    type=values[0],
                    sequence=ref_seq.sequence[int(gene_range[0]):int(gene_range[1])],
                )
                ref_seq.sequences.append(s)

                f = Feature(
                    name=name + '_' + values[0],
                    feature=type[values[0]],
                    start=gene_range[0],
                    end=gene_range[1],
                    strand=strand,
                    attribute='Name=%s;Parent=%s' % (name + '_' + values[0], parent_id),
                    parent_id=parent_id,
                )
                ref_seq.features.append(f)
                s.features.append(f)

        if state and name and values[0] == "3'UTR":
            state = None
            name = None

    db.session.commit()
    return ret
