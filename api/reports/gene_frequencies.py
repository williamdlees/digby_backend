# Gene Frequency plot for AIRR-seq samples
from collections import defaultdict

from werkzeug.exceptions import BadRequest
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file, collate_samples, chunk_list
from app import vdjbase_dbs
from db.vdjbase_model import Gene, GenesDistribution
from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import apply_rep_filter_params, get_multiple_order_file
import pandas as pd

GENE_FREQUENCY_PLOT = 'Gene_Usage.R'
SAMPLE_CHUNKS = 400


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format not in (['pdf', 'html']):
        raise BadRequest('Invalid format requested')

    single_sample_filter = 1 if params['single_sample'] == 'One Selected Sample' else 0
    calc_by_clone = 1 if params['calculate_by'] == 'Number of Clones' else 0
    chain, samples_by_dataset = collate_samples(rep_samples)

    # Format we need to produce is [(gene_name, hetero count, homo count),...]

    genes_frequencies = defaultdict(list)

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].get_session()

        for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id)\
                .filter(Sample.sample_name.in_(sample_chunk))\
                .filter(Sample.sample_group >= single_sample_filter)\
                .all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            i = 0
            sample_list_len = len(sample_list)

            frequencies = session.query(GenesDistribution.sample_id, Gene.name, GenesDistribution.frequency)\
                .join(Gene)\
                .join(Sample)\
                .filter(GenesDistribution.count_by_clones == calc_by_clone)\
                .filter(Gene.name.in_(wanted_genes)) \
                .filter(Sample.sample_name.in_(sample_list)) \
                .all()

            for frequency in frequencies:
                genes_frequencies[frequency[1]].append(round(float(frequency[2]), 2))

    if len(genes_frequencies) == 0:
        raise BadRequest('No frequencies were found for the selected genes')

    labels = ['GENE', 'FREQ']

    genes_frequencies_df = pd.DataFrame(columns=labels)
    #genes_frequencies_df = pd.concat([genes_frequencies_df, pd.DataFrame([{'GENE': gene, 'FREQ': ",".join([str(x) for x in usages])} for gene, usages in genes_frequencies.items()])], axis=0, ignore_index=True)
    genes_frequencies_df = pd.DataFrame([{'GENE': gene, 'FREQ': ",".join([str(x) for x in usages])} for gene, usages in genes_frequencies.items()])

    input_path = make_output_file('tab')
    genes_frequencies_df.to_csv(input_path, sep='\t', index=False)

    output_path = make_output_file(format)
    attachment_filename = '%s_gene_frequency.pdf' % species

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, samples_by_dataset.keys(), [], locus_order=locus_order)

    cmd_line = ["-i", input_path,
                "-o", output_path,
                "-t", 'T' if format == 'html' else 'F',
                "-c", chain,
                "-g", gene_order_file
                ]

    if run_rscript(GENE_FREQUENCY_PLOT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')

