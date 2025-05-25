# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE
#

# Update a non-IMGT reference set for a specific species
import importlib.util

from receptor_utils import simple_bio_seq as simple
from receptor_utils import number_v
from db.genomic_db import Sequence, Gene, AlleleAliasSets
import os


def get_gene_type(label):
    gene = label.split('*')[0]

    if gene[3] == 'J' or gene[3] == 'V':
        return gene[3]
    elif gene[3] == 'D' and '-' in gene and '_' not in gene or ('C' not in gene and ('TRD' in gene or 'TRB' in gene)):
        return 'D'
    else:
        return 'C'


# Add reference sequences


def update_genomic_ref(session, ref_file, dataset):
    if not os.path.isfile(ref_file):
        return f'No reference file {ref_file}'

    refs = simple.read_fasta(ref_file)
    for name, seq in refs.items():
        # determine gene/allele
        if '*' in name:
            gene_name = name.split('*')[0]
        elif 'LJI.Rh_' in name:               # cirelli format
            name = name.replace('LJI.Rh_', '')
            if name[-2] == '.':
                name = name[:-2] + '*' + name[-1]
            name = name.replace('.', '-')
            if '*' not in name:
                name += '*01'
            gene_name = name.split('*')[0]
        else:
            gene_name = name

        if name[3] == 'V' and '.' not in seq:
            print(f'Error in reference set {ref_file}: V-sequence {name} is not gapped')

        gene_id = session.query(Gene.id).filter(Gene.name==gene_name).one_or_none()

        if not gene_id:
            family = ''
            if '-' in gene_name[4:]:
                family = gene_name[4:].split('-')[0]

            gene_type = gene_name[:3] + get_gene_type(gene_name)

            gene_obj = save_gene(session, gene_name, gene_type, family, 0, 0, False)
            gene_id = gene_obj.id
        else:
            gene_id = gene_id[0]

        # check functionality

        functionality = 'Functional'
        notes = ''

        if 'V' in name:
            _, _, notes = number_v.gap_sequence(seq.replace('.', ''), {name: seq}, {name: seq.replace('.', '')})

            if not notes:
                functionality = 'Functional'
            elif 'Stop' in notes:
                functionality = 'Pseudogene'
            else:
                functionality = 'ORF'
        
        s = Sequence(
            name=name,
            imgt_name='',
            type=find_type(name) if dataset != 'IGHC' else 'C-REGION',
            sequence=seq.replace('.', ''),
            novel=False,
            appearances=0,
            deleted=False,
            gapped_sequence=seq,
            functional=functionality,
            notes=notes,
            gene_id=gene_id,
        )
        session.add(s)

    session.commit()


def find_type(name):
    return name[3] + '-REGION'


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


def read_reference(filename):
    records = simple.read_fasta(filename)
    for name in list(records.keys()):
        if '|' in name:             # assume IMGT convention
            rd = name.split('|')[1]
            records[rd] = records[name]
            del records[name]
    return records


# scan subdirs of reference_dir for alias sets
# add an alias set if we find one or more .fasta files in a subdir
def add_alias_sets(reference_dir, session):
    alias_num = 1

    alleles = session.query(Sequence).all()
    ungapped_seqs = {a.sequence: a for a in alleles}

    if not os.path.isdir(reference_dir):
        return

    for subdir in os.listdir(reference_dir):
        subdir_path = os.path.join(reference_dir, subdir)
        if os.path.isdir(subdir_path):
            fasta_files = [f for f in os.listdir(subdir_path) if f.endswith('.fasta')]
            if fasta_files:
                alias_set = AlleleAliasSets(
                    set_name=subdir,
                    alias_number=alias_num
                )
                session.add(alias_set)
                alias_name = f'alias_{alias_num}'
                print(f'Processing alias set {subdir}')
                alias_num += 1

                for fasta_file in fasta_files:
                    recs = read_reference(os.path.join(reference_dir, subdir_path, fasta_file))

                    for allele, sequence in recs.items():
                        if sequence in ungapped_seqs:
                            similar = ungapped_seqs[sequence]

                            # add the allele name to the alias_name column of the similar allele
                            sn = getattr(similar, alias_name)
                            if sn is None or len(sn) == 0:
                                setattr(similar, alias_name, allele)
                            else:
                                setattr(similar, alias_name, f'{sn}|{allele}')

                session.commit()

    print('Alias sets processed')

 