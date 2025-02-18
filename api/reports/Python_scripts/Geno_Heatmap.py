import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import os
import threading
import warnings
import traceback
from dataclasses import dataclass
from typing import Optional
warnings.filterwarnings('ignore')

pdf_lock = threading.Lock()
html_lock = threading.Lock()

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
        print(f"\nGenerating {total} pages...")
    
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

def process_genotype_data(genotypes):
    geno_db = genotypes[['subject', 'gene', 'GENOTYPED_ALLELES', 'k_diff', 'Freq_by_Clone']]
    geno_db.columns = ['subject', 'gene', 'ALLELES', 'K', 'Freq_by_Clone']
    
    numeric_cols = ['K', 'Freq_by_Clone']
    geno_db[numeric_cols] = geno_db[numeric_cols].apply(pd.to_numeric, errors='coerce')
    
    geno_db['ALLELES'] = geno_db['ALLELES'].str.replace('Deletion', 'Del')
    geno_db = geno_db.assign(ALLELES=geno_db['ALLELES'].str.split(',')).explode('ALLELES')
    
    return geno_db

def create_heatmap_data(geno_db, ordered_genes=None):
    heatmap_data = pd.pivot_table(
        geno_db,
        values='Freq_by_Clone',
        index='subject',
        columns='gene',
        aggfunc='mean',
        fill_value=0
    )
    
    if ordered_genes is not None:
        valid_genes = [g for g in ordered_genes if g in heatmap_data.columns]
        heatmap_data = heatmap_data[valid_genes]
    
    return heatmap_data

def split_dataframe(df, chunk_size):
    return [df.iloc[i:i + chunk_size] for i in range(0, len(df), chunk_size)]

def batch_create_annotations(chunk_samples, chunk_genes, chunk_annotations):
    positions = []
    texts = []
    
    sample_to_idx = {sample: idx for idx, sample in enumerate(chunk_samples)}
    gene_to_idx = {gene: idx for idx, gene in enumerate(chunk_genes)}
    
    for _, row in chunk_annotations.iterrows():
        if row['gene'] in gene_to_idx:
            y_pos = sample_to_idx[row['subject']]
            x_pos = gene_to_idx[row['gene']]
            positions.append((x_pos + 0.5, y_pos + 0.5))
            texts.append(row['ALLELES'])
    
    return positions, texts

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
        print(f"\nGenerating {total} pages...")
    
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

def generate_heatmap_html(genotypes, chain='IGH', ordered_genes=None, lk_cutoff=1, file=None, samples_per_page=50, test_mode=True):
    with html_lock:
        try:
            timer = ProcessTimer.start('HTML Generation', test_mode)
            
            timer.log_step("Processing genotypes")
            geno_db = process_genotype_data(genotypes)
            
            timer.log_step("Creating heatmap data")
            heatmap_data = create_heatmap_data(geno_db, ordered_genes)
            
            # Pre-calculate everything possible
            global_min = float(heatmap_data.values.min())
            global_max = float(heatmap_data.values.max())
            
            # Split data into pages
            data_chunks = split_dataframe(heatmap_data, samples_per_page)
            if test_mode:
                data_chunks = data_chunks[:1]
            
            timer.set_pages(len(data_chunks))
            
            # Create HTML files
            base_name, ext = os.path.splitext(file)
            html_files = []
            
            for i, chunk in enumerate(data_chunks):
                page_file = f"{base_name}_page_{i+1}{ext}"
                html_files.append(page_file)
                
                fig = go.Figure(go.Heatmap(
                    z=chunk.values,
                    x=chunk.columns,
                    y=chunk.index,
                    colorscale='Blues',
                    zmin=global_min,
                    zmax=global_max,
                    colorbar=dict(
                        title='Frequency',
                        titleside='right',
                        thickness=20,
                        len=0.75,
                        x=1.02
                    ),
                    hovertemplate='<b>Sample:</b> %{y}<br><b>Gene:</b> %{x}<br><b>Frequency:</b> %{z:.3f}<extra></extra>',
                    showscale=True,
                    xgap=1,
                    ygap=1
                ))
                
                fig.update_layout(
                    title=dict(
                        text=f'{chain} Gene Haplotype Heatmap - Page {i+1}/{len(data_chunks)}',
                        x=0.5,
                        y=0.98,
                        xanchor='center',
                        yanchor='top',
                        font=dict(size=24)
                    ),
                    width=1500,
                    height=min(800 + (len(chunk) * 15), 2000),
                    xaxis=dict(
                        title='Genes',
                        titlefont=dict(size=14),
                        tickangle=45,
                        showgrid=True,
                        gridcolor='rgba(128, 128, 128, 0.2)',
                        gridwidth=1
                    ),
                    yaxis=dict(
                        title='Samples',
                        titlefont=dict(size=14),
                        showgrid=True,
                        gridcolor='rgba(128, 128, 128, 0.2)',
                        gridwidth=1
                    ),
                    font=dict(family='Arial'),
                    margin=dict(l=150, r=100, t=100, b=100),
                    plot_bgcolor='white'
                )
                
                fig.write_html(
                    page_file,
                    include_plotlyjs='cdn',
                    full_html=True,
                    config={'displayModeBar': True, 'displaylogo': False}
                )
                
                timer.log_page(i + 1)
            
            # Create index page
            index_html = '''
            <html>
            <head>
                <title>Heatmap Pages</title>
                <style>
                    body { font-family: Arial, sans-serif; margin: 40px; }
                    h1 { color: #333; }
                    ul { list-style-type: none; padding: 0; }
                    li { margin: 10px 0; }
                    a { color: #0066cc; text-decoration: none; padding: 5px 10px; }
                    a:hover { background-color: #f0f0f0; }
                </style>
            </head>
            <body>
            <h1>Heatmap Pages</h1>
            <ul>
            '''
            
            with open(file, 'w') as f:
                f.write(index_html)
                for i, html_file in enumerate(html_files):
                    page_name = os.path.basename(html_file)
                    f.write(f'<li><a href="{page_name}">Page {i+1}</a></li>\n')
                f.write('</ul></body></html>')
            
            timer.finish()
            return True
            
        except Exception as e:
            print(f"Error: {str(e)}\n{traceback.format_exc()}")
            return False

def generate_heatmap_pdf(genotypes, chain='IGH', ordered_genes=None, lk_cutoff=1, file=None, samples_per_page=10, test_mode=True):
    with pdf_lock:
        try:
            timer = ProcessTimer.start('PDF Generation', test_mode)
            
            timer.log_step("Processing genotypes")
            geno_db = process_genotype_data(genotypes)
            
            timer.log_step("Creating heatmap data")
            heatmap_data = create_heatmap_data(geno_db, ordered_genes)
            
            # Pre-calculate everything possible
            global_min = float(heatmap_data.values.min())
            global_max = float(heatmap_data.values.max())
            filtered_data = geno_db[geno_db['K'] >= lk_cutoff].copy()
            
            # Split data into pages
            data_chunks = split_dataframe(heatmap_data, samples_per_page)
            if test_mode:
                data_chunks = data_chunks[:1]
            
            timer.set_pages(len(data_chunks))
            
            # Generate PDF with optimized settings
            with PdfPages(file) as pdf:
                for i, chunk in enumerate(data_chunks, 1):
                    # Calculate dimensions
                    n_samples = len(chunk)
                    n_genes = len(chunk.columns)
                    width = max(16, n_genes * 0.3)
                    height = max(10, n_samples * 0.4)
                    
                    # Create figure efficiently
                    fig, ax = plt.subplots(figsize=(width, height))
                    
                    # Create optimized heatmap with grid
                    sns.heatmap(
                        chunk,
                        cmap='Blues',
                        vmin=global_min,
                        vmax=global_max,
                        cbar_kws={
                            'label': 'Frequency',
                            'fraction': 0.046,
                            'pad': 0.04
                        },
                        xticklabels=True,
                        yticklabels=True,
                        ax=ax,
                        rasterized=True,
                        linewidths=0.5,
                        linecolor='gray'
                    )
                    
                    # Set white background
                    ax.set_facecolor('white')
                    ax.grid(True, which='major', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
                    
                    ax.set_title(f'{chain} Gene Haplotype Heatmap - Page {i}/{len(data_chunks)}',
                             pad=20, size=14)
                    ax.set_xlabel('Genes', labelpad=10)
                    ax.set_ylabel('Samples', labelpad=10)
                    
                    plt.setp(ax.get_xticklabels(), 
                            rotation=45,
                            ha='right',
                            rotation_mode='anchor',
                            fontsize=max(6, min(8, 300 / n_genes)))
                    
                    plt.setp(ax.get_yticklabels(),
                            fontsize=max(6, min(8, 300 / n_samples)))
                    
                    # Batch process annotations
                    chunk_samples = chunk.index.tolist()
                    chunk_annotations = filtered_data[filtered_data['subject'].isin(chunk_samples)]
                    
                    cell_height = height / n_samples
                    annotation_fontsize = min(6, cell_height * 72 / 4)
                    
                    positions, texts = batch_create_annotations(
                        chunk_samples, 
                        chunk.columns, 
                        chunk_annotations
                    )
                    
                    for (x, y), text in zip(positions, texts):
                        ax.text(x, y, text, ha='center', va='center', fontsize=annotation_fontsize)
                    
                    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
                    pdf.savefig(fig, bbox_inches='tight', dpi=100)
                    plt.close(fig)
                    
                    timer.log_page(i)
            
            timer.finish()
            return True
            
        except Exception as e:
            print(f"Error: {str(e)}\n{traceback.format_exc()}")
            return False