# Genotype report for a single RepSeq sample

from werkzeug.exceptions import BadRequest
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file, find_primer_translations, translate_primer_alleles, translate_primer_genes
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH
from api.reports.report_utils import check_tab_file
import pandas as pd


MULTIPLE_GENOTYPE_SCRIPT = "html_multiple_genotype_hoverText.R"


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) != 1:
        raise BadRequest('This report processes a single repertoire-derived genotype')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    rep_sample = rep_samples[0]
    html = (format == 'html')

    session = vdjbase_dbs[species][rep_sample['dataset']].session
    primer_trans, gene_subs = find_primer_translations(session)
    p = session.query(Sample.genotype).filter(Sample.name == rep_sample['name']).one_or_none()
    p = p[0].replace('samples/', '')
    sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, rep_sample['dataset'], p)

    if not os.path.isfile(sample_path):
        raise BadRequest('Genotype file for sample %s/%s is missing' % (rep_sample['dataset'], rep_sample['name']))

    sample_path = check_tab_file(sample_path)

    # translate pipeline allele names to VDJbase allele names
    genotype = pd.read_csv(sample_path, sep='\t', dtype=str)

    for col in ['alleles', 'GENOTYPED_ALLELES']:
        genotype[col] = [translate_primer_alleles(x, y, primer_trans) for x, y in zip(genotype['gene'], genotype[col])]

    genotype['gene'] = [translate_primer_genes(x, gene_subs) for x in genotype['gene']]
    sample_path = make_output_file('tsv')
    genotype.to_csv(sample_path, sep='\t', index=False)

    report_path = personal_genotype(rep_sample['name'], sample_path, rep_sample['chain'], html)

    if format == 'pdf':
        attachment_filename = '%s_%s_%s_genotype.pdf' % (species, rep_sample['dataset'], rep_sample['name'])
    else:
        attachment_filename = None

    return send_report(report_path, format, attachment_filename)


def personal_genotype(sample_name, genotype_file, chain, html=True):
    output_path = make_output_file('html' if html else 'pdf')
    file_type = 'T' if html else 'F'
    cmd_line = ["-i", genotype_file,
                "-o", output_path,
                "-t", file_type,
                "--samp", sample_name,
                "-c", chain]

    if run_rscript(MULTIPLE_GENOTYPE_SCRIPT, cmd_line) and os.path.getsize(output_path) > 0:
        return output_path
    else:
        raise BadRequest('No output from report')


