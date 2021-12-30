# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE
#

# Update a non-IMGT reference set for a specific species


from receptor_utils import simple_bio_seq as simple
from db.genomic_db import Sequence
import os


# Add reference sequences
def update_genomic_ref(session, ref_file):
    if not os.path.isfile(ref_file):
        return f'No reference file {ref_file}'

    refs = simple.read_fasta(ref_file)
    for name, seq in refs.items():
        # determine gene/allele
        if '*' in name:
            gene = name.split['*'][0]
        elif '.' in name and len(name.split('.')) == 3:     # cirelli format
            gene = name.split('.')[0] + '.' + name.split('.')[1]
        else:
            gene = name

        if 'V' in name and '.' not in seq:
            print(f'Error in reference set {ref_file}: V-sequence {name} is not gapped')

        s = Sequence(
            name=name,
            gene=gene,
            imgt_name='',
            type=find_type(name),
            sequence=seq.replace('.', ''),
            novel=False,
            deleted=False,
            gapped_sequence=seq,
            functional='F',
        )
        session.add(s)

    session.commit()



def find_type(name):
    region_types = {'IGHV': 'V-REGION', 'IGHD': 'D-REGION', 'IGHJ': 'J-REGION'}

    for k,t in region_types.items():
        if k in name:
            return t

    return ''
