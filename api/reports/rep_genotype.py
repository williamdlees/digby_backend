# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report
from api.reports.report_utils import make_output_file, collate_samples, find_primer_translations, translate_primer_alleles, translate_primer_genes

from api.reports.report_utils import trans_df
from app import app, vdjbase_dbs
from db.vdjbase_model import Gene
from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params, get_multiple_order_file
import pandas as pd


MULTIPLE_GENOTYPE_SCRIPT = "html_multiple_genotype_hoverText.R"


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')
    chain, samples_by_dataset = collate_samples(rep_samples)
    genotypes = []

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        primer_trans, gene_subs = find_primer_translations(session)
        sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(samples_by_dataset[dataset])).all()
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        if len(wanted_genes) > 0:
            for (name, genotype, patient_id) in sample_list:
                sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, genotype.replace('samples/', ''))

                if not os.path.isfile(sample_path):
                    continue

                genotype = pd.read_csv(sample_path, sep='\t', dtype=str)
                genotype = trans_df(genotype)

                # translate pipeline allele names to VDJbase allele names
                for col in ['alleles', 'GENOTYPED_ALLELES']:
                    genotype[col] = [translate_primer_alleles(x, y, primer_trans) for x, y in zip(genotype['gene'], genotype[col])]

                genotype['gene'] = [translate_primer_genes(x, gene_subs) for x in genotype['gene']]
                genotype = genotype[genotype.gene.isin(wanted_genes)]

                subject_name = name if len(samples_by_dataset) == 1 else dataset + '_' + name

                if 'subject' not in genotype.columns.values:
                    genotype.insert(0, 'subject', subject_name)
                else:
                    genotype.subject = subject_name

                genotypes.append(genotype)

    if len(genotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    if len(genotypes) > 20:
        raise BadRequest('Please select at most 20 genotypes, or use the Genotype Heatmap report.')

    geno_path = make_output_file('tsv')
    genotypes = pd.concat(genotypes)
    genotypes.to_csv(geno_path, sep='\t')

    if format == 'pdf':
        attachment_filename = '%s_sampled_genotype.pdf' % (species)
    else:
        attachment_filename = None

    if not params['f_pseudo_genes']:
        pseudo = 'F'
    else:
        pseudo = 'T'

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, samples_by_dataset.keys(), locus_order=locus_order)

    output_path = make_output_file('html' if html else 'pdf')

    file_type = 'T' if html else 'F'
    cmd_line = ["-i", geno_path,
                "-o", output_path,
                "-t", file_type,
                "-g", gene_order_file,
                "-c", chain]

    if run_rscript(MULTIPLE_GENOTYPE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')


