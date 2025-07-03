# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest

from api.reports.report_utils import make_output_file, trans_df, collate_samples, chunk_list, find_primer_translations, translate_primer_alleles, translate_primer_genes
from api.reports.reports import run_rscript, send_report
from app import vdjbase_dbs
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype
from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params, get_multiple_order_file
import pandas as pd
from api.reports.Python_scripts.haplotype_heatmap import hapHeatmap

import time
import threading
import traceback
from dataclasses import dataclass
from typing import Optional
from werkzeug.exceptions import BadRequest

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


#HEATMAP_HAPLOTYPE_SCRIPT = "haplotype_heatmap.R"
SAMPLE_CHUNKS = 400


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    
    timer = ProcessTimer.start('Total Report Generation', test_mode=False)
    
    timer.log_step("Starting data collection")
    
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format != 'pdf':
        raise BadRequest('Invalid format requested')

    html = (format == 'html')

    chain, samples_by_dataset = collate_samples(rep_samples)
    haplotypes = pd.DataFrame()

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        primer_trans, gene_subs = find_primer_translations(session)

        haplos = []
        for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]
            haplo_query = session.query(Sample.sample_name, HaplotypesFile.file)\
                .filter(Sample.sample_name.in_(sample_list))\
                .join(SamplesHaplotype, Sample.id == SamplesHaplotype.samples_id)\
                .filter(SamplesHaplotype.haplotypes_file_id == HaplotypesFile.id)\
                .filter(HaplotypesFile.by_gene == params['haplo_gene'])
            haplos.extend(haplo_query.all())

        for name, filename in haplos:
            sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, filename.replace('samples/', ''))

            if not os.path.isfile(sample_path):
                raise BadRequest('Haplotype file %s is missing.' % (sample_path))

            haplotype = pd.read_csv(sample_path, sep='\t', dtype=str)
            haplotype = trans_df(haplotype)
            haplotype['subject'] = name if len(samples_by_dataset) == 1 else dataset + '_' + name

            # translate pipeline allele names to VDJbase allele names

            col_names = list(haplotype.columns.values)
            for i in (2, 3, 4):
                haplotype[col_names[i]] = [translate_primer_alleles(x, y, primer_trans) for x, y in zip(haplotype['gene'], haplotype[col_names[i]])]

            haplotype['gene'] = [translate_primer_genes(x, gene_subs) for x in haplotype['gene']]
            haplotype = haplotype[haplotype.gene.isin(wanted_genes)]

            haplotypes = pd.concat([haplotypes, haplotype], keys=None, ignore_index=True)[haplotype.columns.tolist()]

    timer.log_step("Data collection completed")
    if len(haplotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    #haplo_path = make_output_file('tsv')
    #haplotypes.to_csv(haplo_path, sep='\t', index=False)

    if not params['f_kdiff'] or params['f_kdiff'] == '':
        params['f_kdiff'] = 0

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, samples_by_dataset.keys(), [], locus_order=locus_order)
    output_path = make_output_file('html' if html else 'pdf')
    
    timer.log_step("Starting heatmap generation")
    genes_order = pd.read_csv(gene_order_file, header=None)[0].tolist()
    if hapHeatmap(haplotypes,chain=chain,genes_order=genes_order,file=output_path):
        timer.finish()
        return send_report(output_path, format, f'{species}_haplotype.{format}')
    else:
        raise BadRequest('No output from report')
    
    




