from collections import defaultdict
from werkzeug.exceptions import BadRequest
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file, collate_samples, chunk_list
from app import vdjbase_dbs
from db.vdjbase_model import Gene, GenesDistribution
from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import apply_rep_filter_params, get_multiple_order_file
import pandas as pd
from api.reports.Python_scripts.Gene_Frequencies import create_gene_frequencies_plot

SAMPLE_CHUNKS = 400

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    single_sample_filter = 1 if params['single_sample'] == 'One Selected Sample' else 0
    calc_by_clone = 1 if params['calculate_by'] == 'Number of Clones' else 0
    chain, samples_by_dataset = collate_samples(rep_samples)

    genes_frequencies = defaultdict(list)

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session

        for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id)\
                .filter(Sample.sample_name.in_(sample_chunk))\
                .filter(Sample.sample_group >= single_sample_filter)\
                .all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

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

    # Prepare data for the new plotting function
    genes_frequencies_df = pd.DataFrame([
        {'Gene': gene, 'Frequency': freq}
        for gene, freqs in genes_frequencies.items()
        for freq in freqs
    ])

    input_path = make_output_file('tab')
    genes_frequencies_df.to_csv(input_path, sep='\t', index=False)

    output_path = make_output_file(format)
    attachment_filename = f'{species}_gene_frequency.{format}'

    # Create the plot using the new function
    create_gene_frequencies_plot(
        input_file=input_path,
        output_file=output_path,
        chain=chain,
        format=format
    )

    if os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')