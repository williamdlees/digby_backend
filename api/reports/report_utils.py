#
# Functions to help report processing
#

import os.path
import pandas as pd
from werkzeug.exceptions import BadRequest
import csv
import tempfile
from db.vdjbase_model import Allele

from app import app


""" for the future """
trans_cols = {
    "GENE": 'gene',
    "ALLELES": 'alleles',
    "COUNTS": 'counts',
    "TOTAL": 'total',
    "NOTE": 'note',
    "KH": 'kh',
    "KD": 'kd',
    "KT": 'kt',
    "KQ": 'kq',
    "K_DIFF": 'k_diff',
    "SUBJECT": 'subject',
    "PRIORS_ROW": 'priors_row',
    "PRIORS_COL": 'priors_col',
    "COUNTS1": 'counts1',
    "COUNTS2": 'counts2',
    "COUNTS3": 'counts3',
    "COUNTS4": 'counts4',
    "K1": 'k1',
    "K2": 'k2',
    "K3": 'k3',
    "K4": 'k4',
}

"""
trans_cols = {
    "gene": 'GENE',
    "alleles": 'ALLELES',
    "counts": 'COUNTS',
    "total": 'TOTAL',
    "note": 'NOTE',
    "kh": 'KH',
    "kd": 'KD',
    "kt": 'KT',
    "kq": 'KQ',
    "k_diff": 'K_DIFF',
    "subject": 'SUBJECT',
    "priors_row": 'PRIORS_ROW',
    "priors_col": 'PRIORS_COL',
    "counts1": 'COUNTS1',
    "counts2": 'COUNTS2',
    "counts3": 'COUNTS3',
    "counts4": 'COUNTS4',
    "k1": 'K1',
    "k2": 'K2',
    "k3": 'K3',
    "k4": 'K4',
}
"""

# Check that a tab file exists. If it does, check that the columns are correctly capitalised
# Fix capitalisation if necessary, returning a corrected file


def check_tab_file(filename, dtype=None):
    if not os.path.isfile(filename):
        raise BadRequest('File %s is missing.' % filename)

    if dtype is not None:
        df = pd.read_csv(filename, dtype=dtype, sep='\t')
    else:
        df = pd.read_csv(filename, sep='\t')

    df = trans_df(df)

    filename = make_output_file(os.path.splitext(filename)[1].replace('.', ''))
    df.to_csv(filename, sep='\t', index=True, na_rep='NA', quoting=csv.QUOTE_NONNUMERIC, index_label=False)

    return filename

def trans_df(df):
    renames = list(set(df.columns.values) & set(trans_cols.keys()))

    if len(renames) > 0:
        trans = {x: trans_cols[x] for x in renames}
        df = df.rename(columns=trans)

    return df


# Make a unique file in the output directory
def make_output_file(format):
    with tempfile.NamedTemporaryFile(suffix='.' + format, dir=app.config['OUTPUT_PATH'], delete=False) as fo:
        output_path = fo.name
        fo.close()
    return output_path


# Collate samples from different datasets and determine chain

def collate_samples(rep_samples):
    samples_by_dataset = {}
    chain = None

    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in samples_by_dataset:
            samples_by_dataset[rep_sample['dataset']] = []
            if chain is None:
                chain = rep_sample['chain']
            elif chain != rep_sample['chain']:
                raise BadRequest('This report requires all samples to be selected from the same chain (IGH, IGK, ...')
        samples_by_dataset[rep_sample['dataset']].append(rep_sample['name'])

    return chain, samples_by_dataset

# split a list into chunks of length n

def chunk_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# functions to translate from pipeline allele names (in the tigger/rabhit files) and vdjbase names

def find_primer_translations(session):
    trans = {}
    alleles = session.query(Allele.name, Allele.pipeline_name).all()

    check23 = []

    for name, pipeline_name in alleles:
        if pipeline_name is not None and len(pipeline_name) > 0:
            for pn in pipeline_name.replace(' ', '').split(','):
                trans[pn] = name

    gene_subs = {}

    for pn, gn in trans.items():
        pg, pa = pn.split('*')
        gg, ga = gn.split('*')

        if pg != gg:
            if pg not in gene_subs:
                n = gg.split('-')[1]
                gene_subs[pg] = [gg.split('-')[0], [n]]
            else:
                n = gg.split('-')[1]
                if n not in gene_subs[pg][1]:
                    gene_subs[pg][1].append(n)
            for a in ga.split('_'):
                if '.' in a:
                    n = a.split('.')[0]
                    if n not in gene_subs[pg][1]:
                        gene_subs[pg][1].append(n)

    for k, v in gene_subs.items():
        gene_subs[k][1].sort()
        gene_subs[k][1] = [str(x) for x in gene_subs[k][1]]
        gene_subs[k] = gene_subs[k][0] + '-' + '/'.join(gene_subs[k][1])

    return trans, gene_subs


def translate_primer_alleles(gene, alleles, primer_trans):
    ret = []
    if alleles is not None and isinstance(alleles, str) and len(alleles) > 0:
        for allele in alleles.replace(' ', '').split(','):
            full_name = gene + '*' + allele
            if full_name in primer_trans and primer_trans[full_name] is not None and '*' in primer_trans[full_name]:
                ret.append(primer_trans[full_name].split('*')[1])
            else:
                ret.append(allele)

    return '' if len(ret) == 0 else ','.join(ret)


def translate_primer_genes(gene, gene_subs):
    return gene_subs[gene] if gene in gene_subs else gene
