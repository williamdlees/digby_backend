import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import colorsys
import numpy as np
import traceback

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

def create_gene_frequencies_plot(input, output_file, chain, format='html', downsample=False):
    timer = ProcessTimer.start('Total Report Generation', test_mode=False)
    try:
        timer.log_step("Starting plot generation")
        df = input
        #print("\n=== Starting plot creation ===")
        #print(f"\nInitial data shape: {df.shape}")
        #print("\nFirst few rows of input data:")
        #print(df.head())
        df['TAG'] = df['Gene'].str[3]
        df['FAMILY'] = df['Gene'].str.split('-').str[0]
        
        original_gene_order = df['Gene'].unique()
        available_tags = sorted(df['TAG'].unique())
        #print("\nUnique genes in dataset:", len(original_gene_order))
        #print("Tags present in data:", available_tags)
        
        if downsample:
            #print("\nBefore downsample shape:", df.shape)
            df = downsample_to_10(df)
            #print("After downsample shape:", df.shape)
        
        #print("\nBefore outlier removal shape:", df.shape)
        df = df[df['Frequency'] < 0.99]
        #print("After outlier removal shape:", df.shape)

        # Calculate the maximum y-axis value and round to next 0.05
        max_freq = df['Frequency'].max()
        y_max = get_nice_upper_limit(max_freq)
        #print(f"\nMaximum frequency: {max_freq:.3f}")
        #print(f"Setting y-axis maximum to: {y_max:.2f}")
        timer.log_step("Data processing complete")
        # =============================
        # ======== HTML OUTPUT ========
        # =============================
        if format == 'html':
            if set(available_tags) == {'J', 'V'}:
                # J and V only: stacked
                fig = make_subplots(
                    rows=2, cols=1,
                    subplot_titles=(f"<b>{chain}J Genes</b>", f"<b>{chain}V Genes</b>"),
                    vertical_spacing=0.15,
                    row_heights=[0.4, 0.6]
                )
                tag_to_position = {'J': (1, 1), 'V': (2, 1)}
                fig.update_layout(height=900)
            elif set(available_tags) == {'D', 'V'}:
                # D and V only: stacked
                fig = make_subplots(
                    rows=2, cols=1,
                    subplot_titles=(f"<b>{chain}D Genes</b>", f"<b>{chain}V Genes</b>"),
                    vertical_spacing=0.15,
                    row_heights=[0.4, 0.6]
                )
                tag_to_position = {'D': (1, 1), 'V': (2, 1)}
                fig.update_layout(height=900)
            elif set(available_tags) == {'V'}:
                # Only V: single plot
                fig = make_subplots(
                    rows=1, cols=1,
                    subplot_titles=(f"<b>{chain}V Genes</b>",),
                    row_heights=[1]
                )
                tag_to_position = {'V': (1, 1)}
                fig.update_layout(height=600)
            else:
                # Default: D + J side-by-side, V below
                fig = make_subplots(
                    rows=2, cols=2,
                    specs=[[{}, {}], [{"colspan": 2}, None]],
                    subplot_titles=(f"<b>{chain}D Genes</b>", f"<b>{chain}J Genes</b>", f"<b>{chain}V Genes</b>"),
                    vertical_spacing=0.15,
                    horizontal_spacing=0.1,
                    row_heights=[0.35, 0.65]
                )
                tag_to_position = {'D': (1, 1), 'J': (1, 2), 'V': (2, 1)}
                fig.update_layout(height=1200)
            for tag in ['D', 'J', 'V']:
                if tag not in available_tags:
                    continue
                
                data = df[df['TAG'] == tag]
                row, col = tag_to_position[tag]
                color_map = get_color_map(data['Gene'])
                genes = [gene for gene in original_gene_order if gene in data['Gene'].unique()]
                
                grouped_data = data.groupby('Gene')
                
                for gene, gene_data in grouped_data:
                    #gene_data = data[data['Gene'] == gene]
                    family = gene.split('-')[0]
                    sample_count = len(gene_data)
                    if sample_count < 50:
                        point_custom_data = np.array([[gene, family, sample_count, freq] for freq in gene_data['Frequency']])
                        point_hover_template = (
                            "<b>Gene Information</b><br>" +
                            "Gene: %{customdata[0]}<br>" +
                            "Family: %{customdata[1]}<br>" +
                            "Number of samples: %{customdata[2]}<br>" +
                            "<br><b>Point Data</b><br>" +
                            "Frequency: %{y:.3f}<br><extra></extra>"
                        )
                        fig.add_trace(
                            go.Scatter(
                                x=[gene] * len(gene_data),
                                y=gene_data['Frequency'],
                                mode='markers',
                                marker=dict(color='black', size=6, line=dict(color='white', width=1)),
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
                        distribution_custom_data = np.array([[gene, family, sample_count, stats['p25'], stats['p50'], stats['p75']] for _ in range(len(gene_data))])
                        distribution_hover_template = (
                            "<b>Gene Information</b><br>" +
                            "Name: %{customdata[0]}<br>" +
                            "Family: %{customdata[1]}<br>" +
                            "Number of samples: %{customdata[2]}<br>" +
                            "<br><b>Distribution Statistics</b><br>" +
                            "25%% of samples below: %{customdata[3]:.3f}<br>" +
                            "50%% of samples below: %{customdata[4]:.3f}<br>" +
                            "75%% of samples below: %{customdata[5]:.3f}<br>" +
                            "<br><b>Current Point</b><br>" +
                            "Frequency: %{y:.3f}<br><extra></extra>"
                        )
                        jittered_y = gene_data['Frequency'] + np.random.normal(0, 0.0005, len(gene_data))
                        jittered_y = np.clip(jittered_y, 0, 0.99)
                        fig.add_trace(
                            go.Violin(
                                x=[gene] * len(gene_data),
                                y=jittered_y,
                                line_color='black',
                                fillcolor=color_map[family],
                                opacity=0.4,
                                showlegend=False,
                                points=False,
                                width=0.8,
                                side='both',
                                spanmode='hard',
                                scalemode='width',
                                bandwidth=0.05,
                                hoverinfo='skip'
                            ), row=row, col=col
                        )
                        fig.add_trace(
                            go.Box(
                                x=[gene] * len(gene_data),
                                y=gene_data['Frequency'],
                                boxpoints=False,
                                line=dict(color='black', width=1.5),
                                fillcolor='rgba(255,255,255,0)',
                                width=0.6,
                                showlegend=False,
                                hoverinfo='skip'
                            ), row=row, col=col
                        )
                        fig.add_trace(
                            go.Scatter(
                                x=[gene] * len(gene_data),
                                y=jittered_y,
                                mode='markers',
                                marker=dict(color='rgba(0,0,0,0)', size=6),
                                opacity=1.0,
                                showlegend=False,
                                customdata=distribution_custom_data,
                                hovertemplate=distribution_hover_template,
                                hoveron='points',
                                hoverinfo='all'
                            ), row=row, col=col
                        )
                fig.update_xaxes(
                    tickangle=90,
                    tickfont=dict(size=10, color='#333333', family="Arial, sans-serif"),
                    tickmode='array',
                    tickvals=genes,
                    ticktext=genes,
                    showticklabels=True,
                    row=row, col=col
                )
            fig.update_layout(
                autosize=True,
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
                hoverlabel=dict(bgcolor="white", font_size=14, font_family="Arial")
            )
            fig.update_yaxes(
                title_text="<b>Frequency</b>",
                title_font=dict(size=16, color='#333333', family="Arial, sans-serif", weight='bold'),
                range=[0, y_max],
                tickfont=dict(size=12, color='#333333', family="Arial, sans-serif"),
                dtick=0.1
            )
            fig.update_xaxes(showgrid=True, gridcolor='#CCCCCC')
            fig.update_yaxes(showgrid=True, gridcolor='#CCCCCC')
            fig.write_html(
                output_file,
                full_html=True,
                include_plotlyjs='cdn',
                config={'responsive': True},
                include_mathjax='cdn'
            )
            #print(f"\nHTML file successfully generated: {output_file}")
            return fig
        # =============================
        # ======== PDF OUTPUT =========
        # =============================
        elif format == 'pdf':
            if set(available_tags) == {'J', 'V'}:
                fig = plt.figure(figsize=(20, 12))
                gs = fig.add_gridspec(2, 1, height_ratios=[0.5, 0.5])
                ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
                axes = {'J': ax1, 'V': ax2}
            else:
                fig = plt.figure(figsize=(20, 16))
                gs = fig.add_gridspec(3, 2, height_ratios=[0.05, 0.35, 0.6], width_ratios=[1, 1])
                ax1, ax2, ax3 = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[2, :])
                axes = {'D': ax1, 'J': ax2, 'V': ax3}
            fig.suptitle(f"Gene Frequencies for {chain}", fontsize=24, fontweight='bold', y=0.98)
            tag_data = {tag: df[df['TAG'] == tag] for tag in ['D', 'J', 'V'] if tag in available_tags}
            unique_genes = df['Gene'].unique()
            color_map = get_color_map(unique_genes)
            for tag, ax in axes.items():
                if tag not in tag_data:
                    continue
                data = tag_data[tag]
                grouped_data = data.groupby('Gene')
                small_gene_groups = []  # collect for scatter
                for i, (gene, gene_data) in enumerate(grouped_data):
                    family = gene.split('-')[0]
                    if len(gene_data) < 50:
                        small_gene_groups.append(gene_data)
                    else:
                        sns.violinplot(
                            x='Gene', y='Frequency', data=gene_data, ax=ax,
                            inner=None, cut=0, color=color_map.get(family, 'gray'), alpha=0.4
                        )
                        #sns.boxplot(
                        #    x='Gene', y='Frequency', data=gene_data, ax=ax,
                        #    width=0.2, color='white', showcaps=True,
                        #    boxprops={'zorder': 2, 'facecolor': 'none'},
                        #    showfliers=False, zorder=2
                        #)
                        
                        ax.boxplot(gene_data['Frequency'], positions=[i], widths=0.2, patch_artist=True,
                        boxprops=dict(facecolor=color_map.get(family, 'gray'), alpha=0.4),
                        medianprops=dict(color='black'))
                
                if small_gene_groups:
                    scatter_df = pd.concat(small_gene_groups)
                    ax.scatter(
                        scatter_df['Gene'], scatter_df['Frequency'],
                        color='black', alpha=1.0, s=30, zorder=3,
                        edgecolors='white', linewidth=0.5
                    )
                
                ax.set_title(f"{chain}{tag} Genes", fontsize=16, fontweight='bold')
                ax.set_xlabel('')
                ax.set_ylim(0, y_max)
                ax.set_ylabel('Frequency', fontweight='bold')
                ax.tick_params(axis='x', rotation=90)
                ax.set_facecolor('white')
                ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
                ax.grid(True, linestyle='--', alpha=0.7)
                ax.set_axisbelow(True)
            plt.tight_layout(rect=[0, 0.03, 1, 0.97])
            fig.patch.set_facecolor('white')
            plt.savefig(output_file, format="pdf", bbox_inches='tight')
            #print(f"PDF file successfully generated: {output_file}")
            return None
    except Exception as e:
        print(f"\nAn error occurred: {str(e)}")
        print("Traceback:")
        print(traceback.format_exc())
        return None