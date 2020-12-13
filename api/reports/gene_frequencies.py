# Gene Frequency plot for AIRR-seq samples
from collections import defaultdict

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report, make_output_file
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, HaplotypesFile, SamplesHaplotype, AllelesSample, Gene, Allele, Patient, AllelesPattern, GenesDistribution
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params
from sqlalchemy import func
import pandas as pd

GENE_FREQUENCY_PLOT = 'Gene_Usage.R'



def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format not in ('pdf'):
        raise BadRequest('Invalid format requested')

    single_sample_filter = 1 if params['single_sample'] == 'One Selected Sample' else 0
    calc_by_clone = 1 if params['calculate_by'] == 'Number of Clones' else 0

    samples_by_dataset = {}
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in samples_by_dataset:
            samples_by_dataset[rep_sample['dataset']] = []
        samples_by_dataset[rep_sample['dataset']].append(rep_sample['name'])

    # Format we need to produce is [(gene_name, hetero count, homo count),...]

    chunk = 30
    genes_frequencies = defaultdict(list)

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        sample_list = session.query(Sample.name, Sample.genotype, Sample.patient_id)\
            .filter(Sample.name.in_(samples_by_dataset[dataset]))\
            .filter(Sample.samples_group >= single_sample_filter)\
            .all()
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
        sample_list = [s[0] for s in sample_list]

        i = 0
        sample_list_len = len(sample_list)

        while i < sample_list_len:
            frequencies = session.query(GenesDistribution.sample_id, Gene.name, GenesDistribution.frequency)\
                .join(Gene)\
                .join(Sample)\
                .filter(GenesDistribution.count_by_clones == calc_by_clone)\
                .filter(Gene.name.in_(wanted_genes)) \
                .filter(Sample.name.in_(sample_list[i:min(sample_list_len, i + chunk)])) \
                .all()

            for frequency in frequencies:
                genes_frequencies[frequency[1]].append(round(float(frequency[2]), 2))

            i += chunk


    labels = ['GENE', 'FREQ']
    genes_frequencies_df = pd.DataFrame(columns=labels)
    for gene, usages in genes_frequencies.items():
        genes_frequencies_df = genes_frequencies_df.append({'GENE': gene, 'FREQ':",".join([str(x) for x in usages])}, ignore_index=True)

    input_path = make_output_file('tab')
    genes_frequencies_df.to_csv(input_path, sep='\t', index=False)

    output_path = make_output_file(format)
    attachment_filename = '%s_gene_frequency.pdf' % species

    cmd_line = ["-i", input_path,
                "-o", output_path,
                "-s", SYSDATA,
                "-t", 'T' if format == 'html' else 'F']

    if run_rscript(GENE_FREQUENCY_PLOT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')

