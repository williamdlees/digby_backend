# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE
#

# Read the reference set downloaded from IMGT
# Read human alleles downloaded from IGPDB


from Bio import SeqIO
import yaml
from app import app, db
import sys
from db.feature_db import Sequence, Species
from db.shared import delete_dependencies

# Replace IMGT records in the database
# This deletes all IMGT records previously in place - so all references and samples will need to be updated afterwards
# A better implementation would update the database with changes, without breaking links from existing features

def update_imgt():
    delete_dependencies(None)
#    return ""
    db.session.query(Sequence).filter(Sequence.novel == 0).delete()
    db.session.commit()
    init_imgt_ref()


# Read IMGT reference file and build the reference and codon usage data, using species defined in the config file
def init_imgt_ref():
    with open('track_imgt_config.yaml', 'r') as fc:
        imgt_config = yaml.load(fc, Loader=yaml.FullLoader)

    species = imgt_config['species']

    imgt_reference_genes = None
    try:
        imgt_reference_genes = read_reference(imgt_config['ogre_ref_file'], species)
    except:
        app.logger.error("Error parsing IMGT reference file: %s" % sys.exc_info()[0])

    imgt_gapped_reference_genes = None
    try:
        imgt_gapped_reference_genes = read_reference(imgt_config['gapped_ogre_ref_file'], species)
    except:
        app.logger.error("Error parsing IMGT gapped file: %s" % sys.exc_info()[0])

    for species, genes in imgt_reference_genes.items():
        sp = db.session.query(Species).filter_by(name=species).one_or_none()

        if not sp:
            sp = Species(name=species)
            db.session.add(sp)

        for name, value in genes.items():
            gapped = imgt_gapped_reference_genes[species].get(name, None)

            s = Sequence(
                name=name,
                imgt_name=name,
                type=value[1],
                sequence=str(str(value[0])),
                novel=False,
                gapped_sequence=str(gapped[0]) if gapped else None,
                species=sp,
            )
            db.session.add(s)

    db.session.commit()


def read_reference(filename, species):
    records = {}

    for sp in species.keys():
        records[species[sp]['alias']] = {}

    for rec in SeqIO.parse(filename, 'fasta'):
        rd = rec.description.split('|')
        if rd[2] in species.keys() and (rd[4] in ['V-REGION', 'D-REGION', 'J-REGION']) and (rec.seq is not None):
            records[species[rd[2]]['alias']][rd[1]] = (rec.seq.lower(), rd[4].replace('REGION', 'GENE'))

    return records

def find_family(gene):
    if '-' not in gene:
        return None

    return gene.split('-')[0][4:]

# find the 1-based index of a nucleotide in a gapped reference sequence, given its index in the ungapped sequence
def find_gapped_index(ind_ungapped, species, gene_name):
    gapped_seq = list(imgt_gapped_reference_genes[species][gene_name])

    ind = 0
    found_nt = 0

    while found_nt < ind_ungapped:
        if gapped_seq[ind] != '.':
            found_nt += 1
        ind += 1

    return ind + 1

# Gap a sequence given the closest gapped reference

def gap_sequence(seq, ref):
    i_seq = iter(list(seq))
    i_ref = iter(list(ref))
    ret = ''

    # if reference is partial at the 5' end, pass nucleotides through until we reach the first 'real' reference position
    five_gapped = True

    try:
        while(True):
            r = next(i_ref)
            if r != '.':
                five_gapped = False
                ret += next(i_seq)
            else:
                if five_gapped:
                    ret += next(i_seq)
                else:
                    ret += '.'
    except StopIteration:
        pass

    # if the sequence is longer than the ref, need to add on trailing nucs

    try:
        while(True):
            ret += next(i_seq)
    except StopIteration:
        pass

    return ret
