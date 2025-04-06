# Genotype report for a single sample

import os
from werkzeug.exceptions import BadRequest
from api.reports.genotypes import process_repseq_genotype, process_genomic_genotype
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file
from app import vdjbase_dbs, genomic_dbs
from api.vdjbase.vdjbase import get_order_file


MULTIPLE_GENOTYPE_SCRIPT = "html_multiple_genotype_hoverText.R"


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) + len(genomic_samples) != 1:
        raise BadRequest('This report processes a single genotype')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')

    if len(rep_samples) > 0:
        sample = rep_samples[0]
        session = vdjbase_dbs[species][sample['dataset']].get_session()
        genotype = process_repseq_genotype(sample['sample_name'], [], session, False)
    else:
        sample = genomic_samples[0]
        sample['pcr_target_locus'] = sample['dataset']
        session = genomic_dbs[species][sample['dataset']].get_session()
        genotype = process_genomic_genotype(sample['sample_name'], [], session, True, False)

    if len(genotype) == 0:
        raise BadRequest('Genotype data for sample %s/%s is not available' % (sample['dataset'], sample['sample_name']))

    sample_path = make_output_file('tsv')
    genotype.to_csv(sample_path, sep='\t', index=False)

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_order_file(species, sample['dataset'], locus_order=locus_order)

    report_path = personal_genotype(sample['sample_name'], sample_path, sample['pcr_target_locus'], gene_order_file, html)

    if format == 'pdf':
        attachment_filename = '%s_%s_%s_genotype.pdf' % (species, sample['dataset'], sample['sample_name'])
    else:
        attachment_filename = None

    return send_report(report_path, format, attachment_filename)


def personal_genotype(sample_name, genotype_file, chain, gene_order_file, html=True):
    output_path = make_output_file('html' if html else 'pdf')
    file_type = 'T' if html else 'F'
    cmd_line = ["-i", genotype_file,
                "-o", output_path,
                "-t", file_type,
                "--samp", sample_name,
                "-g", gene_order_file,
                "-c", chain]

    if run_rscript(MULTIPLE_GENOTYPE_SCRIPT, cmd_line) and os.path.getsize(output_path) > 0:
        return output_path
    else:
        raise BadRequest('No output from report')


