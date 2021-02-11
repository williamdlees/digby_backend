# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest

from api.reports.report_utils import trans_df
from api.reports.reports import SYSDATA, run_rscript, send_report, make_output_file
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, Gene
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params
import pandas as pd


HEATMAP_GENOTYPE_SCRIPT = "genotype_heatmap.R"



def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')

    samples_by_dataset = {}
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in samples_by_dataset:
            samples_by_dataset[rep_sample['dataset']] = []
        samples_by_dataset[rep_sample['dataset']].append(rep_sample['name'])

    genotypes = pd.DataFrame()

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        sample_list = session.query(Sample.name, Sample.genotype, Sample.patient_id).filter(Sample.name.in_(samples_by_dataset[dataset])).all()

        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        if len(wanted_genes) > 0:
            for (name, genotype, patient_id) in sample_list:
                sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, genotype.replace('samples/', ''))

                if not os.path.isfile(sample_path):
                    raise BadRequest('Genotype file for sample %s/%s is missing.' % (dataset, name))

                genotype = pd.read_csv(sample_path, sep='\t', dtype=str)
                genotype = trans_df(genotype)
                genotype = genotype[genotype.gene.isin(wanted_genes)]

                subject_name = name if len(samples_by_dataset) == 1 else dataset + '_' + name

                if 'subject' not in genotype.columns.values:
                    genotype.insert(0, 'subject', subject_name)
                else:
                    genotype.subject = subject_name

                genotypes = genotypes.append(genotype)[genotype.columns.tolist()]

    if len(genotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    geno_path = make_output_file('csv')
    genotypes.to_csv(geno_path, sep='\t')

    if format == 'pdf':
        attachment_filename = '%s_genotype.pdf' % species
    else:
        attachment_filename = None

    output_path = make_output_file('html' if html else 'pdf')
    file_type = 'T' if html else 'F'
    cmd_line = ["-i", geno_path,
                "-o", output_path,
                "-s", SYSDATA,
                "-t", file_type,
                "-k", str(params['f_kdiff'])]

    if run_rscript(HEATMAP_GENOTYPE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')


