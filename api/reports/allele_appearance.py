# allele_appearance.py

from werkzeug.exceptions import BadRequest
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file
from app import vdjbase_dbs, genomic_dbs
from db.vdjbase_model import AllelesSample, Gene, Allele
from db.vdjbase_airr_model import Sample, Patient
from db.genomic_db import Sequence as GenomicSequence, SampleSequence as GenomicSampleSequence, Gene as GenomicGene
from db.genomic_airr_model import Sample as GenomicSample
#import os
from api.vdjbase.vdjbase import apply_rep_filter_params
import pandas as pd
from api.reports.Python_scripts.Allelle_Appearance import create_upset_plot
import time
from datetime import datetime
#import plotly.graph_objects as go
#from plotly.subplots import make_subplots
import itertools

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    """
    Run the allele appearance analysis with safe optimization
    """
    start_time = time.time()
    print(f"\nStarting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    gene_matrices = {}
    
    print("Processing repertoire samples...")
    rep_start = time.time()
    
    # Group samples by dataset
    rep_samples_by_dataset = {}
    for rep_sample in rep_samples:
        dataset = rep_sample['dataset']
        if dataset not in rep_samples_by_dataset:
            rep_samples_by_dataset[dataset] = []
        rep_samples_by_dataset[dataset].append(rep_sample['sample_name'])

    # Process each dataset in batch
    for dataset, sample_names in rep_samples_by_dataset.items():
        session = vdjbase_dbs[species][dataset].session
        
        # Get wanted_genes for first sample
        sample_list, wanted_genes = apply_rep_filter_params(params, [(sample_names[0], None, None)], session)
        if not sample_list:
            continue
            
        # Single query for all samples in dataset
        appearances = session.query(
            Patient.patient_name,
            Gene.name,
            Allele.name,
            Gene.locus_order,
            Gene.alpha_order
        ).filter(
            Sample.id == AllelesSample.sample_id,
            Allele.id == AllelesSample.allele_id,
            Gene.id == Allele.gene_id,
            Patient.id == AllelesSample.patient_id,
            Sample.sample_name.in_(sample_names),
            Gene.name.in_(wanted_genes)
        )
        
        if params['novel_alleles'] == 'Exclude':
            appearances = appearances.filter(Allele.novel == 0)
        if params['ambiguous_alleles'] == 'Exclude':
            appearances = appearances.filter(Allele.is_single_allele == 1)
        
        # Process results efficiently in batches
        results = appearances.all()
        for gene, gene_group in itertools.groupby(sorted(results, key=lambda x: x[1]), key=lambda x: x[1]):
            # Extract all patient names and alleles for this gene at once
            patient_allele_pairs = [(p, a.split('*', 1)[1].upper()) for p, _, a, _, _ in gene_group]
            
            if not patient_allele_pairs:
                continue
                
            # Get unique patients and alleles
            unique_patients = list(set(p for p, _ in patient_allele_pairs))
            unique_alleles = list(set(a for _, a in patient_allele_pairs))
            
            # Create DataFrame at once with zeros
            df = pd.DataFrame(
                False,
                index=unique_patients,
                columns=unique_alleles,
                dtype=bool
            )
            
            # Set all True values at once
            for patient, allele in patient_allele_pairs:
                df.loc[patient, allele] = True
            
            gene_matrices[gene] = df

    print(f"Repertoire processing took {time.time() - rep_start:.2f} seconds")

    print("\nProcessing genomic samples...")
    gen_start = time.time()
    
    # Group genomic samples by dataset
    gen_samples_by_dataset = {}
    for gen_sample in genomic_samples:
        dataset = gen_sample['dataset']
        if dataset not in gen_samples_by_dataset:
            gen_samples_by_dataset[dataset] = []
        gen_samples_by_dataset[dataset].append(gen_sample['sample_name'])

    # Process each genomic dataset in batch
    for dataset, sample_names in gen_samples_by_dataset.items():
        session = genomic_dbs[species][dataset].session
        
        # Get wanted_genes for first sample
        sample_list, wanted_genes = apply_rep_filter_params(params, [(sample_names[0],)], session)
        if not sample_list:
            continue
            
        # Single query for all samples in dataset
        appearances = session.query(
            GenomicSample.sample_name,
            GenomicGene.name,
            GenomicSequence.name,
            GenomicGene.locus_order,
            GenomicGene.alpha_order
        ).filter(
            GenomicSample.id == GenomicSampleSequence.sample_id,
            GenomicSequence.id == GenomicSampleSequence.sequence_id,
            GenomicGene.id == GenomicSequence.gene_id,
            GenomicSequence.type.in_(['V-REGION', 'D-REGION', 'J-REGION']),
            GenomicSample.sample_name.in_(sample_names),
            GenomicGene.name.in_(wanted_genes)
        )

        if params['novel_alleles'] == 'Exclude':
            appearances = appearances.filter(GenomicSequence.novel == 0)
        if not params['f_pseudo_genes']:
            appearances = appearances.filter(GenomicSequence.functional == 'Functional')
            
        # Process results
        for sample_name, gene, allele, locus_order, alpha_order in appearances.all():
            allele = allele.split('*', 1)[1].upper()
            
            if gene not in gene_matrices:
                gene_matrices[gene] = pd.DataFrame(dtype=bool)
            
            if allele not in gene_matrices[gene].columns:
                gene_matrices[gene][allele] = False
            
            if sample_name not in gene_matrices[gene].index:
                gene_matrices[gene].loc[sample_name] = False
                
            gene_matrices[gene].loc[sample_name, allele] = True

    print(f"Genomic processing took {time.time() - gen_start:.2f} seconds")

    # Filter out single allele genes and empty matrices
    print("\nPreparing matrices for plotting...")
    plot_start = time.time()
    
    multi_allele_matrices = {
        gene: matrix 
        for gene, matrix in gene_matrices.items() 
        if len(matrix.columns) > 1 and not matrix.empty
    }

    output_path = make_output_file(format)
            
    print(f"Creating plots for {len(multi_allele_matrices)} genes...")
    if format == 'pdf':
        result = create_upset_plot(multi_allele_matrices, output_path)
    if format == 'html':
        result = create_upset_plot(multi_allele_matrices, output_path)
    print(f"Plot creation took {time.time() - plot_start:.2f} seconds")

    end_time = time.time()
    total_time = end_time - start_time
    
    print(f"Total execution time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    
    if result:
        return send_report(output_path, format, f'{species}_allele_appearance.{format}')
    return None

if __name__ == '__main__':
    pass