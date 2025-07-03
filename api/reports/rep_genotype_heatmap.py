import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import os
import time
import threading
import traceback
from dataclasses import dataclass
from typing import Optional
from werkzeug.exceptions import BadRequest

from api.reports.rep_genotype import fake_gene, process_genomic_genotype, process_repseq_genotype
from api.reports.report_utils import collate_samples, chunk_list, collate_gen_samples, make_output_file
from api.reports.reports import send_report
from app import vdjbase_dbs, genomic_dbs
from db.genomic_airr_model import Sample as GenomicSample
from db.vdjbase_airr_model import Sample
from api.vdjbase.vdjbase import apply_rep_filter_params, get_multiple_order_file
from api.reports.Python_scripts.genoHeatmap_html import generate_heatmap_html
from api.reports.Python_scripts.Geno_Heatmap_pdf import generate_heatmap_pdf
from api.reports.genotypes import unfiltered_process_repseq_genotype

@dataclass
class ProcessTimer:
    process_name: str
    test_mode: bool
    process_id: int
    start_time: float
    last_step_time: float
    total_pages: Optional[int] = None
    
    @classmethod
    def start(cls, name: str, test_mode: bool = False):
        start = time.time()
        return cls(
            process_name=name,
            test_mode=test_mode,
            process_id=threading.get_ident(),
            start_time=start,
            last_step_time=start
        )
    
    def log_step(self, step_name: str):
        current = time.time()
        elapsed = current - self.last_step_time
        total = current - self.start_time
        print(f"{step_name}: {elapsed:.2f}s (Total: {total:.2f}s)")
        self.last_step_time = current
    
    def set_pages(self, total: int):
        self.total_pages = total
        print(f"/nGenerating {total} pages...")
    
    def log_page(self, page_num: int):
        if not self.total_pages:
            return
        print(f"Page {page_num}/{self.total_pages} completed in {time.time() - self.last_step_time:.2f}s")
        self.last_step_time = time.time()
    
    def finish(self):
        print(f"\n{'-'*20}")
        mode = 'TEST' if self.test_mode else 'Production'
        total = time.time() - self.start_time
        print(f"{mode} {self.process_name} completed in {total:.2f}s")
        print(f"{'-'*20}\n")

SAMPLE_CHUNKS = 400


from concurrent.futures import ThreadPoolExecutor, as_completed
from db.genomic_db import Gene as GenomicGene, Sequence as GenomicSequence, SampleSequence as GenomicSampleSequence, SampleSequence as GenomicSampleSequence
from db.genomic_airr_model import Sample as GenomicSample
from db.vdjbase_model import Gene, Allele, AllelesSample
from db.vdjbase_airr_model import Sample

def collect_data(timer, species, rep_samples_by_dataset, gen_samples_by_dataset, params, test_mode):
    genotypes = {}
    all_wanted_genes = set()
    fully_haplotyped = 'Only' in params['geno_hap']

    # Initialize sessions
    sessions = {dataset: vdjbase_dbs[species][dataset].session for dataset in rep_samples_by_dataset}
    genomic_sessions = {dataset: genomic_dbs[species][dataset].session for dataset in gen_samples_by_dataset}
    
    def process_genomic_dataset(dataset, sample_names):
        session = genomic_sessions[dataset]
        samples = session.query(GenomicSample).filter(GenomicSample.sample_name.in_(sample_names)).all()
        return dataset, samples
    
    
    from itertools import chain
    
    def process_rep_dataset(dataset, sample_names):
        session = sessions[dataset]
        
        if len(sample_names) <= SAMPLE_CHUNKS:  # Run a single query if it's small
            results = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id) \
                            .filter(Sample.sample_name.in_(sample_names)).all()
            return dataset, results
        
        # If sample_names is large, process in chunks
        chunks = [sample_names[i:i + SAMPLE_CHUNKS] for i in range(0, len(sample_names), SAMPLE_CHUNKS)]
        
        results_iter = chain.from_iterable(
            session.query(Sample.sample_name, Sample.genotype, Sample.patient_id)
                .filter(Sample.sample_name.in_(chunk)).all()
            for chunk in chunks
        )
        
        return dataset, list(results_iter)  # Convert iterator to a list only at the end

    
    datasets = {}
    for dataset, names in rep_samples_by_dataset.items():
        datasets[dataset] = process_rep_dataset(dataset, names) 
    
    for dataset in datasets:
        chain_type,sample_list = datasets[dataset]
        session = sessions[dataset]
        
        if sample_list:
            timer.log_step("Starting apply_rep_filter_params...")

            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            if test_mode and wanted_genes:
                wanted_genes = list(wanted_genes)[:10]
                
            timer.log_step("Finished apply_rep_filter_params...")
            
            if wanted_genes:
                all_wanted_genes.update(wanted_genes)
                timer.log_step("Starting process_repseq_genotype...")
                if(len(sample_list) < 200):
                    batch_genotypes = {name: process_repseq_genotype(name, all_wanted_genes, session, False)
                                        for name, _, _ in sample_list}
                else:
                    timer.log_step("Querying Full Table...")
                    full_table = unfiltered_process_repseq_genotype(wanted_genes, session, False)
                    timer.log_step("Full Table Here...")
                    sample_ids = set(s[0] for s in sample_list)  # set is faster for lookup
                    filtered = full_table[full_table['subject'].isin(sample_ids)]
                    batch_genotypes = dict(tuple(filtered.groupby('subject')))
                
                genotypes.update(batch_genotypes)

            timer.log_step("Finished process_repseq_genotype...")
            
    timer.log_step("Process rep_seq data in parallel completed")
    
    # Process genomic data in parallel
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_dataset = {executor.submit(process_genomic_dataset, dataset, names): dataset
                             for dataset, names in gen_samples_by_dataset.items()}

        for future in as_completed(future_to_dataset):
            dataset = future_to_dataset[future]
            samples = future.result()[1]
            if not samples:
                continue

            session = genomic_sessions[dataset]
            sample_list = [(s.sample_name, s.annotation_path, s.sample_name) for s in samples]
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

            if test_mode and wanted_genes:
                wanted_genes = list(wanted_genes)[:10]

            if wanted_genes:
                all_wanted_genes.update(wanted_genes)
                filtered_samples = [s[0] for s in sample_list]
                batch_genotypes = {name: process_genomic_genotype(name, all_wanted_genes, session,
                                                                  not params['f_pseudo_genes'], fully_haplotyped)
                                   for name in filtered_samples}
                genotypes.update(batch_genotypes)

    timer.log_step("Process genomic data in parallel completed")
    
    timer.log_step("Data collection completed")
    return genotypes, all_wanted_genes


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params, test_mode=False):
    timer = ProcessTimer.start('Total Report Generation', test_mode)
    
    html = (format == 'html')
    chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)

    if not chain:
        chain = g_chain

    if test_mode:
        for dataset in rep_samples_by_dataset:
            rep_samples_by_dataset[dataset] = rep_samples_by_dataset[dataset][:5]
        for dataset in gen_samples_by_dataset:
            gen_samples_by_dataset[dataset] = gen_samples_by_dataset[dataset][:5]

    timer.log_step("Starting data collection")
    
    genotypes, all_wanted_genes = collect_data(timer, species, rep_samples_by_dataset, 
                                             gen_samples_by_dataset, params, test_mode)

    if not genotypes:
        raise BadRequest('No records matching the filter criteria were found.')

    timer.log_step("Processing genotypes")
    for subject_name, genotype in genotypes.items():
        contained_genes = genotype['gene'].tolist() if len(genotype) > 0 else []
        missing_genes = all_wanted_genes - set(contained_genes)
        if missing_genes:
            fake_genes = [fake_gene({'alleles': ['Unk'], 'count': [], 'fc': [], 'fs': [], 
                                   'total_count': 0, 'kdiff': 0}, gene, subject_name)
                         for gene in missing_genes]
            genotypes[subject_name] = pd.concat([genotype, pd.DataFrame(fake_genes)], 
                                              ignore_index=True)
    
    for genotype in genotypes.values():
        genotype.sort_values(by=['gene'], inplace=True)

    #geno_path = make_output_file('tsv')
    combined_genotypes = pd.concat(genotypes.values())
    #combined_genotypes.to_csv(geno_path, sep='\t')

    timer.log_step("Starting plot generation")
    output_path = make_output_file('html' if html else 'pdf')
    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, rep_samples_by_dataset.keys(), 
                                            gen_samples_by_dataset.keys(), locus_order=locus_order)

    if create_report(combined_genotypes, output_path, html, chain, gene_order_file, params['f_kdiff']):
        timer.finish()
        return send_report(output_path, format, f'{species}_genotype.{format}')
    else:
        raise BadRequest('No output from report')
    
    timer.log_step("plot generation completed")

def create_report(genotypes, output_path, html, chain, gene_order_file, f_kdiff):
    output_dir = os.path.dirname(output_path)
    output_name = os.path.basename(output_path)
    
    gene_order = None
    if gene_order_file and os.path.exists(gene_order_file):
        gene_order = pd.read_csv(gene_order_file, sep='\t', header=None).iloc[:, 0].tolist()

    try:
        if html:
            return generate_heatmap_html(
                genotypes,
                chain=chain,
                gene_sort=gene_order,
                file=os.path.join(output_dir, output_name)
            )
        else:
            return generate_heatmap_pdf(
                genotypes,
                chain=chain,
                gene_sort=gene_order,
                file=os.path.join(output_dir, output_name)
            )
    except Exception as e:
        raise BadRequest(f'Error generating report: {str(e)}')