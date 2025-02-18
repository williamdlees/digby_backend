import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import colorsys
import numpy as np
import traceback

# Optional: Set some general matplotlib parameters for better plotting
plt.style.use('default')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

def generate_colors(n):
    HSV_tuples = [(x * 1.0 / n, 0.5, 0.9) for x in range(n)]
    return list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))

def get_color_map(genes):
    families = sorted(set([gene.split('-')[0] for gene in genes]))
    colors = generate_colors(len(families))
    return {family: f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' for family, (r, g, b) in zip(families, colors)}

def get_distribution_stats(data):
    """Calculate distribution statistics for a gene's frequency data"""
    percentiles = np.percentile(data, [25, 50, 75])
    return {
        'p25': round(percentiles[0], 3),
        'p50': round(percentiles[1], 3),
        'p75': round(percentiles[2], 3)
    }

def downsample_to_10(df):
    def sample_gene(group):
        if len(group) > 10:
            return group.sample(n=10, random_state=42)
        return group
    return df.groupby('Gene').apply(sample_gene).reset_index(drop=True)

def get_nice_upper_limit(max_val):
    """Get a nice upper limit for the y-axis that's a multiple of 0.05"""
    return np.ceil(max_val * 20) / 20  # This will round up to the next 0.05

def create_gene_frequencies_plot(input_file, output_file, chain, format='html', downsample=False):
    try:
        print("\n=== Starting plot creation ===")
        df = pd.read_csv(input_file, sep='\t')
        print(f"\nInitial data shape: {df.shape}")
        print("\nFirst few rows of input data:")
        print(df.head())
        
        df['TAG'] = df['Gene'].str[3]
        df['FAMILY'] = df['Gene'].str.split('-').str[0]
        
        original_gene_order = df['Gene'].unique()
        print("\nUnique genes in dataset:", len(original_gene_order))
        
        if downsample:
            print("\nBefore downsample shape:", df.shape)
            df = downsample_to_10(df)
            print("After downsample shape:", df.shape)
        
        print("\nBefore outlier removal shape:", df.shape)
        df = df[df['Frequency'] < 0.99]
        print("After outlier removal shape:", df.shape)

        # Calculate the maximum y-axis value and round to next 0.05
        max_freq = df['Frequency'].max()
        y_max = get_nice_upper_limit(max_freq)
        print(f"\nMaximum frequency: {max_freq:.3f}")
        print(f"Setting y-axis maximum to: {y_max:.2f}")

        if format == 'html':
            fig = make_subplots(
                rows=3, cols=2,
                subplot_titles=(f"<b>{chain}D Genes</b>", f"<b>{chain}J Genes</b>", None, None, None, None),
                vertical_spacing=0.1,
                horizontal_spacing=0.1,
                row_heights=[0.3, 0.1, 0.6]
            )

            for i, tag in enumerate(['D', 'J', 'V']):
                print(f"\nProcessing {tag} genes...")
                data = df[df['TAG'] == tag]
                print(f"Number of {tag} genes: {len(data['Gene'].unique())}")
                
                if tag in ['D', 'J']:
                    row, col = (1, i+1)
                else:
                    row, col = (3, 1)
                
                color_map = get_color_map(data['Gene'])
                genes = [gene for gene in original_gene_order if gene in data['Gene'].unique()]
                
                for gene in genes:
                    gene_data = data[data['Gene'] == gene]
                    family = gene.split('-')[0]
                    sample_count = len(gene_data)
                    is_low_sample = sample_count < 50
                    print(f"\nProcessing gene {gene}: {sample_count} samples")
                    
                    if is_low_sample:
                        point_custom_data = np.array([[
                            gene,
                            family,
                            sample_count,
                            freq
                        ] for freq in gene_data['Frequency']])
                        
                        point_hover_template = (
                            "<b>Gene Information</b><br>" +
                            "Gene: %{customdata[0]}<br>" +
                            "Family: %{customdata[1]}<br>" +
                            "Number of samples: %{customdata[2]}<br>" +
                            "<br>" +
                            "<b>Point Data</b><br>" +
                            "Frequency: %{y:.3f}<br>" +
                            "<extra></extra>"
                        )
                        
                        fig.add_trace(
                            go.Scatter(
                                x=[gene] * len(gene_data),
                                y=gene_data['Frequency'],
                                mode='markers',
                                marker=dict(
                                    color='black',
                                    size=6,
                                    line=dict(color='white', width=1),
                                    symbol='circle'
                                ),
                                opacity=1.0,
                                showlegend=False,
                                customdata=point_custom_data,
                                hovertemplate=point_hover_template,
                                hoveron='points'
                            ),
                            row=row, col=col
                        )
                    else:
                        stats = get_distribution_stats(gene_data['Frequency'])
                        
                        distribution_custom_data = np.array([[
                            gene,
                            family,
                            sample_count,
                            stats['p25'],
                            stats['p50'],
                            stats['p75']
                        ] for _ in range(len(gene_data))])

                        distribution_hover_template = (
                            "<b>Gene Information</b><br>" +
                            "Name: %{customdata[0]}<br>" +
                            "Family: %{customdata[1]}<br>" +
                            "Number of samples: %{customdata[2]}<br>" +
                            "<br>" +
                            "<b>Distribution Statistics</b><br>" +
                            "25%% of samples below: %{customdata[3]:.3f}<br>" +
                            "50%% of samples below: %{customdata[4]:.3f}<br>" +
                            "75%% of samples below: %{customdata[5]:.3f}<br>" +
                            "<br>" +
                            "<b>Current Point</b><br>" +
                            "Frequency: %{y:.3f}<br>" +
                            "<extra></extra>"
                        )

                        jittered_y = gene_data['Frequency'] + np.random.normal(0, 0.0005, len(gene_data))
                        jittered_y = np.clip(jittered_y, 0, 0.99)
                        
                        # Add violin plot
                        fig.add_trace(
                            go.Violin(
                                x=[gene] * len(gene_data),
                                y=jittered_y,
                                name=gene,
                                box_visible=False,
                                line_color='black',
                                meanline_visible=False,
                                fillcolor=color_map[family],
                                opacity=0.4,  # Reduced opacity for better boxplot visibility
                                showlegend=False,
                                points=False,
                                width=0.8,
                                side='both',
                                spanmode='hard',
                                scalemode='width',
                                bandwidth=0.05,
                                hoverinfo='skip',
                                line_width=1
                            ), 
                            row=row, col=col
                        )

                        # Add box plot
                        fig.add_trace(
                            go.Box(
                                x=[gene] * len(gene_data),
                                y=gene_data['Frequency'],
                                name=gene,
                                boxpoints=False,  # Hide individual points
                                line=dict(color='black', width=1.5),
                                fillcolor='rgba(255,255,255,0)',  # Transparent fill
                                width=0.6,  # Slightly narrower than violin plot
                                showlegend=False,
                                hoverinfo='skip'
                            ),
                            row=row, col=col
                        )

                        # Add transparent scatter points with hover template
                        fig.add_trace(
                            go.Scatter(
                                x=[gene] * len(gene_data),
                                y=jittered_y,
                                mode='markers',
                                marker=dict(
                                    color='rgba(0,0,0,0)',
                                    size=6
                                ),
                                opacity=1.0,
                                showlegend=False,
                                customdata=distribution_custom_data,
                                hovertemplate=distribution_hover_template,
                                hoveron='points',
                                hoverinfo='all'
                            ),
                            row=row, col=col
                        )
                
                fig.update_xaxes(
                    tickangle=90,
                    tickfont=dict(size=10, color='#333333', family="Arial, sans-serif"),
                    tickmode='array',
                    tickvals=genes,
                    ticktext=genes,
                    title_font=dict(size=16, color='#333333', family="Arial, sans-serif", weight='bold'),
                    showticklabels=True,
                    row=row, col=col
                )

            fig.update_layout(
                autosize=True,
                height=1200,
                title=dict(
                    text=f"<b>Gene Frequencies for {chain} (10 Samples per Gene)</b>" if downsample else f"<b>Gene Frequencies for {chain}</b>",
                    font=dict(size=28, color='#333333'),
                    x=0.5,
                    xanchor='center'
                ),
                margin=dict(b=50, t=100, l=50, r=50),
                plot_bgcolor='white',
                paper_bgcolor='white',
                font=dict(family="Arial, sans-serif"),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=14,
                    font_family="Arial"
                )
            )

            # Update subplot titles and add centered IGHV title
            for i, annotation in enumerate(fig['layout']['annotations']):
                if i < 2:
                    annotation.update(
                        font=dict(size=18, color='#333333', family="Arial, sans-serif", weight='bold')
                    )
                else:
                    fig['layout']['annotations'].pop(i)

            fig.add_annotation(
                x=0.5,
                y=0.62,
                xref="paper",
                yref="paper",
                text=f"<b>{chain}V Genes</b>",
                showarrow=False,
                font=dict(size=20, color='#333333', family="Arial, sans-serif", weight='bold'),
                xanchor='center',
                yanchor='bottom'
            )

            # Update y-axes with rounded limits and nice tick values
            fig.update_yaxes(
                title_text="<b>Frequency</b>",
                title_font=dict(size=16, color='#333333', family="Arial, sans-serif", weight='bold'),
                range=[0, y_max],
                tickfont=dict(size=12, color='#333333', family="Arial, sans-serif"),
                dtick=0.1  # This creates ticks at 0, 0.1, 0.2, etc.
            )

            fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='#CCCCCC')
            fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='#CCCCCC')

            fig.update_layout(
                xaxis=dict(domain=[0, 0.45]),
                xaxis2=dict(domain=[0.55, 1]),
                xaxis5=dict(domain=[0, 1]),
                yaxis=dict(domain=[0.7, 1]),
                yaxis2=dict(domain=[0.7, 1]),
                yaxis5=dict(domain=[0, 0.6])
            )

            fig.write_html(
                output_file,
                full_html=True,
                include_plotlyjs='cdn',
                config={'responsive': True},
                include_mathjax='cdn'
            )

            print(f"\nHTML file successfully generated: {output_file}")
            return fig
        
        elif format == 'pdf':
            with PdfPages(output_file) as pdf:
                fig = plt.figure(figsize=(20, 16))
                gs = fig.add_gridspec(3, 2, height_ratios=[0.05, 0.35, 0.6], width_ratios=[1, 1])
                
                fig.suptitle(f"Gene Frequencies for {chain}", fontsize=24, fontweight='bold', y=0.98)
                
                ax1 = fig.add_subplot(gs[1, 0])
                ax2 = fig.add_subplot(gs[1, 1])
                ax3 = fig.add_subplot(gs[2, :])
                
                for tag, ax in zip(['D', 'J', 'V'], [ax1, ax2, ax3]):
                    data = df[df['TAG'] == tag]
                    color_map = get_color_map(data['Gene'])
                    
                    for gene in data['Gene'].unique():
                        gene_data = data[data['Gene'] == gene]
                        family = gene.split('-')[0]
                        
                        if len(gene_data) < 50:
                            x = np.full(len(gene_data), gene_data['Gene'].iloc[0])
                            y = gene_data['Frequency']
                            ax.scatter(x=x, y=y, color='black', alpha=1.0, s=30, 
                                    zorder=3, edgecolors='white', linewidth=0.5)
                        else:
                            # Add violin plot
                            sns.violinplot(x='Gene', y='Frequency', data=gene_data, ax=ax, 
                                        inner=None, cut=0, color=color_map[family], alpha=0.4)
                            # Add box plot overlay
                            sns.boxplot(x='Gene', y='Frequency', data=gene_data, ax=ax,
                                    width=0.2, color='white', showcaps=True, 
                                    boxprops={'zorder': 2, 'facecolor': 'none'},
                                    showfliers=False, zorder=2)

                    ax.set_title(f"{chain}{tag} Genes", fontsize=16, fontweight='bold')
                    ax.set_xlabel('')
                    ax.set_ylim(0, y_max)  # Use rounded upper limit
                    ax.set_ylabel('Frequency', fontweight='bold')
                    ax.tick_params(axis='x', rotation=90)
                    ax.set_facecolor('white')
                    
                    # Set major ticks at 0.1 intervals
                    ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
                    
                    ax.grid(True, linestyle='--', alpha=0.7)
                    ax.set_axisbelow(True)

                plt.tight_layout(rect=[0, 0.03, 1, 0.97])
                fig.patch.set_facecolor('white')
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

                print(f"PDF file successfully generated: {output_file}")
            return None

    except Exception as e:
        print(f"\nAn error occurred: {str(e)}")
        print("Traceback:")
        print(traceback.format_exc())
        return None