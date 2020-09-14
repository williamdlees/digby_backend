# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest

from api.reports.report_utils import trans_df
from api.reports.reports import SYSDATA, run_rscript, send_report, make_output_file
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, HaplotypesFile, SamplesHaplotype
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params
import pandas as pd

HEATMAP_HAPLOTYPE_SCRIPT = "haplotype_heatmap.R"


def run(format, species, genomic_samples, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format != 'pdf':
        raise BadRequest('Invalid format requested')

    html = (format == 'html')

    samples_by_dataset = {}
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in samples_by_dataset:
            samples_by_dataset[rep_sample['dataset']] = []
        samples_by_dataset[rep_sample['dataset']].append(rep_sample['name'])

    haplotypes = pd.DataFrame()

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        sample_list = session.query(Sample.name, Sample.genotype, Sample.patient_id).filter(Sample.name.in_(samples_by_dataset[dataset])).all()
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
        sample_list = [s[0] for s in sample_list]
        haplos = session.query(Sample.name, HaplotypesFile.file)\
            .filter(Sample.name.in_(sample_list))\
            .join(SamplesHaplotype)\
            .filter(Sample.id == SamplesHaplotype.samples_id)\
            .filter(SamplesHaplotype.haplotypes_files_id == HaplotypesFile.id)\
            .filter(HaplotypesFile.by_gene == params['haplo_gene']).all()

        for name, filename in haplos:
            sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, filename.replace('samples/', ''))

            if not os.path.isfile(sample_path):
                raise BadRequest('Haplotype file %s is missing.' % (sample_path))

            haplotype = pd.read_csv(sample_path, sep='\t', dtype=str)
            haplotype = trans_df(haplotype)
            haplotype = haplotype[haplotype.GENE.isin(wanted_genes)]
            haplotype['SUBJECT'] = name if len(samples_by_dataset) == 1 else dataset + '_' + name
            haplotypes = pd.concat([haplotypes, haplotype], keys=None, ignore_index=True)[haplotype.columns.tolist()]

    if len(haplotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    haplo_path = make_output_file('csv')
    haplotypes.to_csv(haplo_path, sep='\t', index=False)
    attachment_filename = '%s_haplotype_heatmap.pdf' % species

    if not params['f_kdiff'] or params['f_kdiff'] == '':
        params['f_kdiff'] = 0

    output_path = make_output_file('html' if html else 'pdf')
    cmd_line = ["-i", haplo_path,
                "-o", output_path,
                "-s", SYSDATA,
                "-k", str(params['f_kdiff'])]

    if run_rscript(HEATMAP_HAPLOTYPE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')


