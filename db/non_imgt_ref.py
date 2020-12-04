# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE
#

# Update a non-IMGT reference set for a specific species


from Bio import SeqIO
from app import app, db
from db.feature_db import Sequence, Species
import os


# Add non IMGT reference sequences
def update_non_imgt_ref(species, ref_file):
    if not os.path.isfile(ref_file):
        return 'No reference file %s for species %s' % (ref_file, species)

    sp = db.session.query(Species).filter_by(name=species).one_or_none()

    if not sp:
        sp = Species(name=species)
        db.session.add(sp)

    with open(ref_file, 'r') as fi:
        refs = SeqIO.parse(fi, 'fasta')

        for ref in refs:
            s = Sequence(
                name=ref.name,
                imgt_name='',
                type=find_type(ref.name),
                sequence=str(ref.seq).lower(),
                novel=False,
                deleted=False,
                gapped_sequence='',
                species=sp,
                functional='F',
            )
            db.session.add(s)

        db.session.commit()

    return 'Processed non-IMGT reference set for species %s' % species


def find_type(name):
    region_types = {'IGHV': 'V-REGION', 'IGHD': 'D-REGION', 'IGHJ': 'J-REGION'}

    for k,t in region_types.items():
        if k in name:
            return t

    return ''
