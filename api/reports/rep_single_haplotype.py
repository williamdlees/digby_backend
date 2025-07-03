# Haplotype report for a single RepSeq sample

from werkzeug.exceptions import BadRequest

from api.reports.reports import SYSDATA, run_rscript, send_report
from api.reports.report_utils import make_output_file, check_tab_file, find_primer_translations, translate_primer_alleles, translate_primer_genes
import pandas as pd
from app import vdjbase_dbs
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype
from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, get_order_file
from api.reports.Python_scripts.haplotype_plot import generate_haplotype_plot


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) != 1:
        raise BadRequest('This report processes a single repertoire-derived haplotype')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    rep_sample = rep_samples[0]
    html = (format == 'html')
    print(rep_sample)

    session = vdjbase_dbs[species][rep_sample['dataset']].session
    primer_trans, gene_subs = find_primer_translations(session)

    p = session.query(HaplotypesFile.file)\
        .join(SamplesHaplotype)\
        .join(Sample)\
        .filter(Sample.sample_name == rep_sample['sample_name'])\
        .filter(HaplotypesFile.by_gene_s == params['haplo_gene']).one_or_none()
    print(p)
    p = p[0].replace('samples/','')
    sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, rep_sample['dataset'], p)
    print(p)
    if not os.path.isfile(sample_path):
        raise BadRequest('Genotype file for sample %s/%s is missing' % (rep_sample['dataset'], rep_sample['sample_name']))

    sample_path = check_tab_file(sample_path)

    # translate pipeline allele names to VDJbase allele names
    haplotype = pd.read_csv(sample_path, sep='\t', dtype=str)

    col_names = list(haplotype.columns.values)
    for i in (2, 3, 4):
        haplotype[col_names[i]] = [translate_primer_alleles(x, y, primer_trans) for x, y in zip(haplotype['gene'], haplotype[col_names[i]])]

    haplotype['gene'] = [translate_primer_genes(x, gene_subs) for x in haplotype['gene']]
    #sample_path = make_output_file('tsv')
    #haplotype.to_csv(sample_path, sep='\t', index=False)

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_order_file(species, rep_sample['dataset'], locus_order=locus_order)

    report_path = personal_haplotype(rep_sample['sample_name'], haplotype, rep_sample['pcr_target_locus'], gene_order_file, html)

    if format == 'pdf':
        attachment_filename = '%s_%s_%s_%s_haplotype.pdf' % (species, rep_sample['dataset'], rep_sample['sample_name'], params['haplo_gene'])
    else:
        attachment_filename = None

    return send_report(report_path, format, attachment_filename)

def personal_haplotype(sample_name, haplotype, chain, gene_order_file, html=True):
    output_path = make_output_file('html' if html else 'pdf')
    output_dir = os.path.dirname(output_path)
    output_name = os.path.basename(output_path)
    
    gene_order = None
    if gene_order_file and os.path.exists(gene_order_file):
        gene_order = pd.read_csv(gene_order_file, sep='\t', header=None).iloc[:, 0].tolist()
    try:
        return generate_haplotype_plot(
            sample_name = sample_name,
            hap_table=haplotype,
            chain=chain,
            genes_order=gene_order,
            file=os.path.join(output_dir, output_name),
            html = html
        )
    except Exception as e:
        raise BadRequest(f'Error generating report: {str(e)}')



