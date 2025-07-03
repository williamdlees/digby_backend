#import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
#import plotly.io as pio
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#from matplotlib.backends.backend_pdf import PdfPages
#import io
#from PIL import Image

def create_html_plot(df, chain):
    """Generates the heterozygous plot dynamically using Plotly with properly formatted hover text."""

    gene_segment = df.copy()

    # Normalize counts to create ratios
    gene_segment['Total'] = gene_segment['HM'] + gene_segment['HT']
    gene_segment['HT_Ratio'] = gene_segment['HT'] / gene_segment['Total']
    gene_segment['HM_Ratio'] = gene_segment['HM'] / gene_segment['Total']

    # Extract gene types dynamically based on the first four characters
    gene_segment['GeneType'] = gene_segment['GENE'].str[:4]
    gene_types = sorted(gene_segment['GeneType'].unique())  # Sort for consistency

    # Define colors
    colors = {'Homozygous': '#A7C7E7', 'Heterozygous': '#F6A6A6'}

    # Create subplots dynamically based on the number of unique gene types
    fig = make_subplots(
        rows=1, cols=len(gene_types),
        subplot_titles=[f'<b>{gene_type} Genes</b>' for gene_type in gene_types],
        horizontal_spacing=0.05,
        column_widths=[1/len(gene_types)] * len(gene_types)
    )

    for i, gene_type in enumerate(gene_types, start=1):
        data = gene_segment[gene_segment['GeneType'] == gene_type]

        # Corrected hover text formatting without HTML styles
        ht_hover_text = [
            f"Gene: {gene}<br>"
            f"Ratio: {ht_ratio:.2f}<br>"
            f"{ht} from {total} subjects<br>"
            f"──────────────<br>"
            f"Heterozygous"
            for gene, ht_ratio, ht, total in zip(data['GENE'], data['HT_Ratio'], data['HT'], data['Total'])
        ]

        hm_hover_text = [
            f"Gene: {gene}<br>"
            f"Ratio: {hm_ratio:.2f}<br>"
            f"{hm} from {total} subjects<br>"
            f"──────────────<br>"
            f"Homozygous"
            for gene, hm_ratio, hm, total in zip(data['GENE'], data['HM_Ratio'], data['HM'], data['Total'])
        ]

        fig.add_trace(
            go.Bar(
                y=data['GENE'],
                x=data['HT_Ratio'],
                name='Heterozygous',
                marker_color=colors['Heterozygous'],
                orientation='h',
                showlegend=i == 1,
                customdata=ht_hover_text,
                hovertemplate="%{customdata}<extra></extra>",
                width=0.7
            ),
            row=1, col=i
        )

        fig.add_trace(
            go.Bar(
                y=data['GENE'],
                x=data['HM_Ratio'],
                name='Homozygous',
                marker_color=colors['Homozygous'],
                orientation='h',
                showlegend=i == 1,
                base=data['HT_Ratio'],
                customdata=hm_hover_text,
                hovertemplate="%{customdata}<extra></extra>",
                width=0.7
            ),
            row=1, col=i
        )

        fig.update_yaxes(title_text='<b>Gene</b>' if i == 1 else None, row=1, col=i)
        fig.update_xaxes(title_text='<b>Ratio</b>', row=1, col=i, range=[0, 1])

    for i in range(1, len(gene_types) + 1):
        fig.update_yaxes(tickfont=dict(size=8), row=1, col=i) 
    
    fig.update_layout(
        title_text=f'<b>Heterozygosity Distribution for Chain {chain}</b>',
        title_x=0.5,
        height=1200,
        width=1500,
        barmode='stack',
        plot_bgcolor='white',
        legend=dict(
            x=1.05,  # Move legend to the right
            y=0.5,
            traceorder="normal",
            font=dict(size=12)
        )
    )

    return fig

def create_pdf_plot(df, chain):

    """Generates the heterozygous plot with dynamic gene types extracted from gene names."""

    # Load the dataset
    gene_segment = df

    # Normalize counts to create ratios
    gene_segment['Total'] = gene_segment['HM'] + gene_segment['HT']
    gene_segment['HT_Ratio'] = gene_segment['HT'] / gene_segment['Total']
    gene_segment['HM_Ratio'] = gene_segment['HM'] / gene_segment['Total']

    # Extract gene types dynamically based on the 4th character of 'GENE'
    gene_segment['GeneType'] = gene_segment['GENE'].str[:4]
    gene_types = sorted(gene_segment['GeneType'].unique())  # Sort for consistency

    # Define colors
    colors = {'Homozygous': '#A7C7E7', 'Heterozygous': '#F6A6A6'}

    # Determine max number of genes for dynamic figure sizing
    max_genes = max(gene_segment[gene_segment['GeneType'] == gt].shape[0] for gt in gene_types)

    # Dynamically adjust figure size
    fig_width = 6 * len(gene_types)  # 6 inches per subplot
    fig_height = max(5, min(15, max_genes * 0.5))  # Scale with gene count, min 5, max 15

    # Create figure with subplots
    fig, axes = plt.subplots(1, len(gene_types), figsize=(fig_width, fig_height), sharey=False)
    fig.suptitle(f'Heterozygosity Distribution for Chain {chain} \n', fontsize=16, fontweight='bold')

    # Ensure axes is iterable even for a single subplot
    if len(gene_types) == 1:
        axes = [axes]

    for i, (gene_type, ax) in enumerate(zip(gene_types, axes)):
        data = gene_segment[gene_segment['GeneType'] == gene_type].copy()

        # Sort data to ensure a consistent order
        data = data.sort_values(by='HT_Ratio', ascending=True)

        # Get values for plotting
        y_labels = data['GENE']
        HT_Ratio = data['HT_Ratio']
        HM_Ratio = data['HM_Ratio']

        # Plot heterozygous (HT) bar
        ax.barh(y_labels, HT_Ratio, color=colors['Heterozygous'], label='Heterozygous')

        # Plot homozygous (HM) bar stacked on HT
        ax.barh(y_labels, HM_Ratio, color=colors['Homozygous'], label='Homozygous', left=HT_Ratio)

        # Set axis labels
        ax.set_xlim(0, 1)
        ax.set_xlabel('Ratio', fontsize=10, fontweight='bold')
        ax.set_title(f'{gene_type} Genes', fontsize=12, fontweight='bold')
        ax.grid(axis='x', linestyle='--', alpha=0.7)

    # Move legend to the side (right)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center', bbox_to_anchor=(1, 0.5), fontsize=12)

    plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout for the legend
    return plt
    
def create_heterozygous_plot_html(df, output_html, chain):
    """Creates the heterozygous plot and saves it as an HTML file."""
    fig = create_html_plot(df, chain)
    fig.write_html(output_html)

def create_heterozygous_plot_pdf(df, output_pdf, chain):
    plt = create_pdf_plot(df, chain)
    #fig.write_image(output_pdf)
    plt.savefig(output_pdf, format="pdf", bbox_inches='tight')