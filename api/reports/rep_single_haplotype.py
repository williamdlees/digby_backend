# Haplotype report for a single RepSeq sample

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report, make_output_file
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, HaplotypesFile, SamplesHaplotype
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH


PERSONAL_HAPLOTYPE_SCRIPT = 'Haplotype_plot.R'


def run(format, species, genomic_samples, rep_samples, params):
    if len(rep_samples) != 1:
        raise BadRequest('This report processes a single repertoire-derived haplotype')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    rep_sample = rep_samples[0]
    html = (format == 'html')

    session = vdjbase_dbs[species][rep_sample['dataset']].session
    p = session.query(HaplotypesFile.file)\
        .join(SamplesHaplotype)\
        .join(Sample)\
        .filter(Sample.name == rep_sample['name'])\
        .filter(HaplotypesFile.by_gene_s == params['haplo_gene']).one_or_none()

    p = p[0].replace('samples/','')
    sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, rep_sample['dataset'], p)
    report_path = personal_haplotype(rep_sample['name'], sample_path, html)

    if format == 'pdf':
        attachment_filename = '%s_%s_%s_%s_haplotype.pdf' % (species, rep_sample['dataset'], rep_sample['name'], params['haplo_gene'])
    else:
        attachment_filename = None

    return send_report(report_path, format, attachment_filename)


def personal_haplotype(sample_name, haplotype_file, html=True):
    output_path = make_output_file('html' if html else 'pdf')
    file_type = 'T' if html else 'F'
    cmd_line = ["-i", haplotype_file,
                "-o", output_path,
                "-s", SYSDATA,
                "-t", file_type,
                "--samp", sample_name]

    if run_rscript(PERSONAL_HAPLOTYPE_SCRIPT, cmd_line) and os.path.getsize(output_path) > 0:
        return output_path
    else:
        raise BadRequest('No output from report')


