# Genotype report for a single sample

import os
from werkzeug.exceptions import BadRequest
from api.reports.genotypes import process_repseq_genotype, process_genomic_genotype
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file
from app import vdjbase_dbs, genomic_dbs
from api.vdjbase.vdjbase import get_order_file
import pandas as pd
from api.reports.Python_scripts.genotype_plot import generate_genotype_plot


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) + len(genomic_samples) != 1:
        raise BadRequest('This report processes a single genotype')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')

    if len(rep_samples) > 0:
        sample = rep_samples[0]
        session = vdjbase_dbs[species][sample['dataset']].session
        genotype = process_repseq_genotype(sample['sample_name'], [], session, False)
    else:
        sample = genomic_samples[0]
        sample['pcr_target_locus'] = sample['dataset']
        session = genomic_dbs[species][sample['dataset']].session
        genotype = process_genomic_genotype(sample['sample_name'], [], session, True, False)

    if len(genotype) == 0:
        raise BadRequest('Genotype data for sample %s/%s is not available' % (sample['dataset'], sample['sample_name']))

    #sample_path = make_output_file('tsv')
    #genotype.to_csv(sample_path, sep='\t', index=False)

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_order_file(species, sample['dataset'], locus_order=locus_order)

    report_path = personal_genotype(sample['sample_name'], genotype, sample['pcr_target_locus'], gene_order_file, html)

    if format == 'pdf':
        attachment_filename = '%s_%s_%s_genotype.pdf' % (species, sample['dataset'], sample['sample_name'])
    else:
        attachment_filename = None

    return send_report(report_path, format, attachment_filename)

def personal_genotype(sample_name, genotypes, chain, gene_order_file, html=True):
    output_path = make_output_file('html' if html else 'pdf')
    output_dir = os.path.dirname(output_path)
    output_name = os.path.basename(output_path)
    
    gene_order = None
    if gene_order_file and os.path.exists(gene_order_file):
        gene_order = pd.read_csv(gene_order_file, sep='\t', header=None).iloc[:, 0].tolist()
    try:
        return generate_genotype_plot(
            genotypes,
            chain=chain,
            gene_sort=gene_order,
            file=os.path.join(output_dir, output_name),
            html = html
        )
    except Exception as e:
        raise BadRequest(f'Error generating report: {str(e)}')



