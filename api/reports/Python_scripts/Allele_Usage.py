import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots

def create_allele_usage_plot(gene_segment, chain="IGH"):
    # Validate input
    if not isinstance(gene_segment, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame")
    if not all(col in gene_segment.columns for col in ['GENE', 'COUNT']):
        raise ValueError("Input DataFrame must have 'GENE' and 'COUNT' columns")
    if chain not in ["IGH", "IGK", "IGL", "TRB", "TRA"]:
        raise ValueError("Invalid chain value")

    # Define colors
    pastel_blue = '#BAE1FF'
    background_color = '#F8F8F8'  # Light gray background
    grid_color = '#E0E0E0'  # Slightly darker gray for grid lines

    # Determine gene segment (V, D, or J)
    chain_prefix = chain[0]
    nth = 4 if gene_segment['GENE'].str.startswith(chain_prefix).any() else 1
    gene_segment['SEGMENT'] = gene_segment['GENE'].str[nth - 1]
    
    # Get unique segments
    unique_segments = gene_segment['SEGMENT'].unique()
    
    if len(unique_segments) == 0:
        raise ValueError("No segments found in the data")

    # Create subplot for each segment
    fig = make_subplots(rows=1, cols=len(unique_segments), shared_xaxes=False, 
                        subplot_titles=[f"<b>{chain}{seg}</b>" for seg in unique_segments])

    # Plot data for each segment
    for i, segment in enumerate(unique_segments):
        segment_data = gene_segment[gene_segment['SEGMENT'] == segment]
        if not segment_data.empty:
            trace = go.Bar(
                y=segment_data['GENE'],
                x=segment_data['COUNT'],
                orientation='h',
                marker=dict(color=pastel_blue),
                name=f"{chain}{segment}"
            )
            fig.add_trace(trace, row=1, col=i+1)
        else:
            print(f"No data available for segment {segment}")

    # Update layout with improved aesthetics, centered headline, and enhanced grid
    fig.update_layout(
        height=1300,
        showlegend=False,
        title=dict(
            text=f"<b>Allele Usage for {chain}</b>",
            font=dict(family="Arial, sans-serif", size=24, color="#000000"),
            x=0.5,
            xanchor='center'
        ),
        barmode='stack',
        font=dict(family="Arial, sans-serif", size=12, color="#000000"),
        plot_bgcolor=background_color,
        paper_bgcolor='white'
    )

    # Update x-axes
    fig.update_xaxes(
        title_text="<b>Count</b>", 
        title_font=dict(size=14), 
        tickfont=dict(size=12),
        showgrid=True,
        gridcolor=grid_color,
        gridwidth=1
    )

    # Update y-axes
    fig.update_yaxes(
        title_text="<b>Gene</b>", 
        title_font=dict(size=14), 
        tickfont=dict(size=12),
        showgrid=False
    )

    return fig