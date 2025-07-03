# Allele usage plot for AIRR-seq samples

# Import necessary libraries and modules
from werkzeug.exceptions import BadRequest
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file, collate_samples, chunk_list, collate_gen_samples

from app import vdjbase_dbs, genomic_dbs, app
from db.vdjbase_model import AllelesSample, Gene, Allele, AllelesPattern
from db.genomic_db import Sequence as GenomicSequence, SampleSequence as GenomicSampleSequence, Gene as GenomicGene
from db.genomic_airr_model import Sample as GenomicSample

from db.vdjbase_airr_model import Patient, Sample
import os
from api.vdjbase.vdjbase import apply_rep_filter_params
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import time
from api.reports.Python_scripts.Allele_Usage import create_allele_usage_plot, create_allele_usage_plot_pdf


# Define constant for sample chunk size
SAMPLE_CHUNKS = 400

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    
    start_time = time.time()
    
    # Validate input format
    if format not in ["html", "pdf"]:
        raise ValueError("Invalid format. Choose 'html' or 'pdf'.")
    
    # Process input parameters
    kdiff = float(params['f_kdiff']) if 'f_kdiff' in params and params['f_kdiff'] != '' else 0
    chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)
    # Validate chain consistency
    if (chain and g_chain) and chain != g_chain:
        raise BadRequest('This report requires all samples to be selected from the same chain (IGH, IGK, ...')

    if not chain:
        chain = g_chain

    # Validate ambiguous alleles processing
    if len(rep_samples_by_dataset) + len(gen_samples_by_dataset) > 1 and params['ambiguous_alleles'] != 'Exclude':
        raise BadRequest('Ambiguous alleles cannot be processed across multiple datasets')

    # Initialize dictionary to store gene allele counts
    gene_allele_counts = {}

    # Process repertoire samples
    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        allele_recs = []

        for sample_chunk in chunk_list(rep_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            # Query and filter samples
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            # Build query for alleles
            query = session.query(Gene.name, Allele.name, Gene.type) \
                .join(Allele) \
                .join(AllelesSample) \
                .join(Sample) \
                .join(Patient, Patient.id == Sample.patient_id) \
                .filter(Gene.name.in_(wanted_genes)) \
                .filter(Allele.name.notlike('%Del%')) \
                .filter(Allele.name.notlike('%OR%')) \
                .filter(Sample.sample_name.in_(sample_list)) \
                .filter(AllelesSample.kdiff >= kdiff)

            # Apply sorting
            if 'sort_order' in params and params['sort_order'] == 'Locus':
                query = query.order_by(Gene.locus_order, Patient.id, Allele.name)
            else:
                query = query.order_by(Gene.alpha_order, Patient.id, Allele.name)

            # Apply filters for novel and ambiguous alleles
            if params['novel_alleles'] == 'Exclude':
                query = query.filter(Allele.novel == 0)

            if params['ambiguous_alleles'] == 'Exclude':
                query = query.filter(Allele.is_single_allele == 1)

            allele_recs.extend(query.all())

        # Process allele records
        i = 0
        while i < len(allele_recs):
            (gene_name, allele_name, gene_type) = allele_recs[i]
            gene_allele_names = []

            while i < len(allele_recs):
                if allele_recs[i][0] != gene_name:
                    break

                allele_name = allele_recs[i][1]
                gene_allele_names.append(allele_name)
                i += 1

            gene_allele_names = set(gene_allele_names)

            # Handle ambiguous alleles
            if(params['ambiguous_alleles'] != 'Exclude'):
                patterns = session.query(AllelesPattern.pattern_id)\
                    .filter(AllelesPattern.allele_in_p_id.in_(gene_allele_names))\
                    .filter(AllelesPattern.pattern_id.in_(gene_allele_names))\
                    .all()

                if patterns is not None and len(patterns) > 0:
                    patterns = set([pattern[0] for pattern in patterns])
                    gene_allele_names = gene_allele_names - patterns

            # Update gene allele counts
            if gene_name not in gene_allele_counts:
                gene_allele_counts[gene_name] = gene_allele_names
            else:
                gene_allele_counts[gene_name] |= gene_allele_names

    # Process genomic samples
    for dataset in gen_samples_by_dataset.keys():
        session = genomic_dbs[species][dataset].session
        allele_recs = []

        for sample_chunk in chunk_list(gen_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            # Query and filter samples
            sample_list = session.query(GenomicSample.sample_name).filter(GenomicSample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            # Build query for alleles
            query = session.query(GenomicGene.name, GenomicSequence.name, GenomicGene.type) \
                .filter(GenomicSample.id == GenomicSampleSequence.sample_id) \
                .filter(GenomicSequence.id == GenomicSampleSequence.sequence_id) \
                .filter(GenomicGene.id == GenomicSequence.gene_id)\
                .filter(GenomicSequence.type.in_(['V-REGION', 'D-REGION', 'J-REGION'])) \
                .filter(GenomicSample.sample_name.in_(sample_list))\
                .filter(GenomicGene.name.in_(wanted_genes))

            # Apply sorting
            if 'sort_order' in params and params['sort_order'] == 'Locus':
                query = query.order_by(GenomicGene.locus_order, GenomicSample.sample_name, GenomicSequence.name)
            else:
                query = query.order_by(GenomicGene.alpha_order, GenomicSample.sample_name, GenomicSequence.name)

            # Apply filters for novel alleles and pseudo genes
            if params['novel_alleles'] == 'Exclude':
                query = query.filter(GenomicSequence.novel == 0)

            if not params['f_pseudo_genes']:
                query = query.filter(GenomicSequence.functional == 'Functional')

            allele_recs.extend(query.all())

        # Process allele records
        i = 0
        while i < len(allele_recs):
            (gene_name, allele_name, gene_type) = allele_recs[i]
            gene_allele_ids = []

            while i < len(allele_recs):
                if allele_recs[i][0] != gene_name:
                    break

                allele_name = allele_recs[i][1]
                gene_allele_ids.append(allele_name)
                i += 1

            gene_allele_ids = set(gene_allele_ids)

            # Update gene allele counts
            if gene_name not in gene_allele_counts:
                gene_allele_counts[gene_name] = gene_allele_ids
            else:
                gene_allele_counts[gene_name] |= gene_allele_ids

    # Prepare data for plotting
    listed_allele_count = []
    for gene, alleles in gene_allele_counts.items():
        listed_allele_count.append((gene, len(alleles)))

    # Create DataFrame
    labels = ['GENE', 'COUNT']
    df = pd.DataFrame(listed_allele_count, columns=labels)
    
    # Generate and save report
    output_path = make_output_file(format)
    attachment_filename = f'{species}_allele_usage.{format}'
    generate_report(format,df, output_path, chain)

    # Log execution time
    end_time = time.time()
    execution_time = end_time - start_time
    app.logger.info(f"The report creation took {execution_time:.6f} seconds")

    # Return report or raise exception if no output
    if os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format,attachment_filename)
    else:
        raise BadRequest('No output from report')

def allele_usage_bar_html(gene_segment, chain="IGH"):
    # Validate chain
    if chain not in ["IGH", "IGK", "IGL", "TRB", "TRA"]:
        raise ValueError("Invalid chain value")

    # Determine gene segment (V, D, or J)
    chain_prefix = chain[0]
    nth = 4 if gene_segment['GENE'].str.startswith(chain_prefix).any() else 1
    gene_segment['SEGMENT'] = gene_segment['GENE'].apply(lambda x: x[nth - 1])
    
    # Get unique segments
    unique_segments = gene_segment['SEGMENT'].unique()
    
    # Create subplot for each segment
    fig = make_subplots(rows=1, cols=len(unique_segments), shared_xaxes=False, subplot_titles=[f"{chain}{seg}" for seg in unique_segments])

    # Plot data for each segment
    for i, segment in enumerate(unique_segments):
        segment_data = gene_segment[gene_segment['SEGMENT'] == segment]
        segment_data = segment_data.sort_values(by='GENE', ascending=False)
        if not segment_data.empty:
            trace = go.Bar(
                y=segment_data['GENE'],
                x=segment_data['COUNT'],
                orientation='h',
                marker=dict(color='peru'),
                name=f"{chain}{segment}"
            )
            fig.add_trace(trace, row=1, col=i+1)
        else:
            print(f"No data available for segment {segment}")

    # Update layout
    fig.update_layout(height=1300, showlegend=False, title_text=f"Allele Usage for {chain}", barmode='stack')
    return fig

def generate_report(format, df, output_file, chain="IGH"):
    # Generate plot 
    if format == "html":
        fig = create_allele_usage_plot(df, chain)
        fig.write_html(output_file)
    elif format == "pdf":
        try:
            plt=create_allele_usage_plot_pdf(df, chain)
            plt.savefig(output_file, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Error generating PDF: {e}")
