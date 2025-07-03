#!/usr/bin/env python

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate genotype heatmaps')
    
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='Input TSV file containing genotype data')
    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help='Output file name for the generated graph')
    parser.add_argument('-t', '--is_html', type=str, required=True,
                        help='Type of output file: "T" for HTML, "F" for PDF')
    parser.add_argument('-k', '--Kdiff', type=float, required=True,
                        help='The minimal kdiff value')
    parser.add_argument('-c', '--chain', type=str, default='IGH',
                        help='Chain type: IGH, IGK, IGL, TRB, TRA')
    parser.add_argument('-g', '--gene_order_file', type=str, default=None,
                        help='TSV file listing genes in the desired order')
    
    return parser.parse_args()


def read_genotype_data(file_path):
    """Read genotype data from TSV file."""
    return pd.read_csv(file_path, sep='\t')


def read_gene_order(file_path):
    """Read gene order from TSV file if provided."""
    if file_path:
        gene_order_df = pd.read_csv(file_path, sep='\t', header=None)
        return gene_order_df[0].tolist()
    return None


def generate_heatmap_pdf(genotypes, chain, ordered_genes, lk_cutoff, output_file):
    """Generate a PDF heatmap using matplotlib/seaborn."""
    # Prepare data for heatmap
    # Assuming genotypes dataframe has columns: gene, sample, value, etc.
    
    # Filter based on lk_cutoff if needed
    if 'lk_diff' in genotypes.columns:
        genotypes = genotypes[genotypes['lk_diff'] >= lk_cutoff]
    
    # Pivot data for heatmap format
    pivot_data = genotypes.pivot_table(
        index='sample', 
        columns='gene', 
        values='value',  # Adjust this column name based on your actual data
        aggfunc='first'
    )
    
    # Reorder genes if ordered_genes is provided
    if ordered_genes:
        # Get intersection of available genes and ordered genes
        available_genes = set(pivot_data.columns)
        ordered_available = [g for g in ordered_genes if g in available_genes]
        
        # Reorder columns
        pivot_data = pivot_data[ordered_available]
    
    # Create figure with appropriate size
    num_genes = len(pivot_data.columns)
    num_samples = len(pivot_data.index)
    plt.figure(figsize=(max(8, num_genes * 0.24 + 1.5), max(6, num_samples * 0.2 + 1.5)))
    
    # Generate heatmap
    sns.heatmap(pivot_data, cmap='viridis', linewidths=0.5)
    
    plt.title(f'Genotype Heatmap - {chain}')
    plt.tight_layout()
    
    # Save to file
    plt.savefig(output_file)
    plt.close()


def generate_heatmap_html(genotypes, chain, ordered_genes, lk_cutoff, output_file):
    """Generate an interactive HTML heatmap using Plotly."""
    # Filter based on lk_cutoff if needed
    if 'lk_diff' in genotypes.columns:
        genotypes = genotypes[genotypes['lk_diff'] >= lk_cutoff]
    
    # Pivot data for heatmap format
    pivot_data = genotypes.pivot_table(
        index='sample', 
        columns='gene', 
        values='value',  # Adjust based on your actual data
        aggfunc='first'
    )
    
    # Reorder genes if ordered_genes is provided
    if ordered_genes:
        # Get intersection of available genes and ordered genes
        available_genes = set(pivot_data.columns)
        ordered_available = [g for g in ordered_genes if g in available_genes]
        
        # Reorder columns
        pivot_data = pivot_data[ordered_available]
    
    # Create figure
    fig = px.imshow(
        pivot_data,
        labels=dict(x="Gene", y="Sample", color="Value"),
        title=f'Genotype Heatmap - {chain}',
        color_continuous_scale='Viridis',
        aspect="auto"
    )
    
    fig.update_layout(
        height=max(600, len(pivot_data.index) * 20),
        width=max(800, len(pivot_data.columns) * 24 + 150),
        xaxis_nticks=len(pivot_data.columns),
        yaxis_nticks=len(pivot_data.index)
    )
    
    # Save as HTML
    fig.write_html(output_file)


def main():
    """Main function to execute the heatmap generation."""
    args = parse_arguments()
    
    # Validate required arguments
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")
    
    # Read genotype data
    genotypes = read_genotype_data(args.input_file)
    
    # Read gene order if provided
    gene_order = read_gene_order(args.gene_order_file) if args.gene_order_file else None
    
    # Normalize output path
    output_dir = os.path.normpath(os.path.dirname(args.output_file))
    output_file = os.path.join(output_dir, os.path.basename(args.output_file))
    
    # Convert is_html string to boolean
    is_html = args.is_html.upper() == 'T'
    
    # Generate heatmap based on output type
    if is_html:
        generate_heatmap_html(genotypes, args.chain, gene_order, args.Kdiff, output_file)
    else:
        generate_heatmap_pdf(genotypes, args.chain, gene_order, args.Kdiff, output_file)
    
    print(f"Heatmap saved to: {output_file}")


if __name__ == "__main__":
    main()