# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE
#

# Update a non-IMGT reference set for a specific species
import importlib

from receptor_utils import simple_bio_seq as simple
from db.genomic_db import Sequence, Gene
import os


# Add reference sequences


def update_genomic_ref(session, ref_file):
    if not os.path.isfile(ref_file):
        return f'No reference file {ref_file}'

    refs = simple.read_fasta(ref_file)
    for name, seq in refs.items():
        # determine gene/allele
        if '*' in name:
            gene_name = name.split('*')[0]
        elif 'LJI.Rh_' in name:               # cirelli format
            if name[-2:-1] == '.' and name[-1].isalpha():
                gene_name = name[:-2]
            else:
                gene_name = name
        else:
            gene_name = name

        if 'V' in name and '.' not in seq:
            print(f'Error in reference set {ref_file}: V-sequence {name} is not gapped')

        gene_id = session.query(Gene.id).filter(Gene.name==gene_name).one_or_none()

        if not gene_id:
            family = ''
            if '-' in gene_name[4:]:
                family = gene_name[4:].split('-')[0]


            gene_type = gene_name.replace('LJI.Rh_', '')[:4]

            gene_obj = save_gene(session, gene_name, gene_type, family, 0, 0, False)
            gene_id = gene_obj.id
        else:
            gene_id = gene_id[0]

        s = Sequence(
            name=name,
            imgt_name='',
            type=find_type(name),
            sequence=seq.replace('.', ''),
            novel=False,
            appearances=0,
            deleted=False,
            gapped_sequence=seq,
            functional='Functional',
            notes='',
            gene_id=gene_id,
        )
        session.add(s)

    session.commit()



def find_type(name):
    region_types = {'IGHV': 'V-REGION', 'IGHD': 'D-REGION', 'IGHJ': 'J-REGION'}

    for k,t in region_types.items():
        if k in name:
            return t

    return ''


def read_gene_order(session, dataset_dir):
    if os.path.isfile(os.path.join(dataset_dir, 'gene_order.py')):
        spec = importlib.util.spec_from_file_location("gene_order", os.path.join(dataset_dir, 'gene_order.py'))
        gene_order = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gene_order)
    else:
        print('gene_order.py not found - skipped')
        return

    alpha_order = 0
    for gene in gene_order.ALPHA_ORDER:
        locus_order = gene_order.LOCUS_ORDER.index(gene) if gene in gene_order.LOCUS_ORDER else 999
        pseudo = gene in gene_order.PSEUDO_GENES

        family = ''
        if '-' in gene[4:]:
            family = gene[4:].split('-')[0]

        gene_type = gene[:4]

        save_gene(session, gene, gene_type, family, locus_order, alpha_order, pseudo)
        alpha_order += 1

    session.commit()


def save_gene(session, name, type, family, locus_order, alpha_order, pseudo_gene):
    g = Gene(
        name=name,
        type=type,
        family=family,
        locus_order=locus_order,
        alpha_order=alpha_order,
        pseudo_gene=pseudo_gene,
    )
    session.add(g)
    session.flush()
    return g
