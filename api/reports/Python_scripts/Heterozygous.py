import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def create_heterozygous_plot(input_file, output_file, chain):
    # Load the data
    gene_segment = pd.read_csv(input_file, sep='\t')

    # Normalize counts to create ratios
    gene_segment['Total'] = gene_segment['HM'] + gene_segment['HT']
    gene_segment['HT_Ratio'] = gene_segment['HT'] / gene_segment['Total']
    gene_segment['HM_Ratio'] = gene_segment['HM'] / gene_segment['Total']

    # Calculate the number of genes in each segment
    gene_counts = gene_segment['GENE'].apply(lambda gene: gene[:4]).value_counts()

    subplot_height = 800

    # Create subplots, one for each gene type, side by side
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('<b>IGHV Genes</b>', '<b>IGHD Genes</b>', '<b>IGHJ Genes</b>'),
        horizontal_spacing=0.05,
        column_widths=[0.4, 0.3, 0.3]
    )

    gene_types = ['IGHV', 'IGHD', 'IGHJ']
    colors = {'Homozygous': '#A7C7E7', 'Heterozygous': '#F6A6A6'}

    for i, gene_type in enumerate(gene_types, start=1):
        # Filter data for the current gene type
        data = gene_segment[gene_segment['GENE'].str.startswith(gene_type)]
        

        # Add bar traces for Heterozygous ratio
        fig.add_trace(
            go.Bar(
                y=data['GENE'],
                x=data['HT_Ratio'],
                name='Heterozygous',
                marker_color=colors['Heterozygous'],
                orientation='h',
                showlegend=i == 1,
                hovertemplate='<b>Gene:</b> %{y}<br><b>HT Count:</b> %{customdata[0]:.0f} out of %{customdata[1]:.0f}<br><b>Ratio:</b> %{x:.2%}<extra></extra>',
                customdata=data[['HT', 'Total']].values,
                width=0.75
            ),
            row=1, col=i
        )

        # Add bar traces for Homozygous ratio
        fig.add_trace(
            go.Bar(
                y=data['GENE'],
                x=data['HM_Ratio'],
                name='Homozygous',
                marker_color=colors['Homozygous'],
                orientation='h',
                showlegend=i == 1,
                hovertemplate='<b>Gene:</b> %{y}<br><b>HM Count:</b> %{customdata[0]:.0f} out of %{customdata[1]:.0f}<br><b>Ratio:</b> %{customdata[2]:.2%}<extra></extra>',
                base=data['HT_Ratio'],
                customdata=data[['HM', 'Total', 'HM_Ratio']].values,
                width=0.75
            ),
            row=1, col=i
        )

        # Update y-axis for better gene name visibility
        fig.update_yaxes(
            title_text='<b>Gene</b>' if i == 1 else None,
            title_font=dict(size=15, color='black', family='Arial Black'),
            row=1, col=i,
            tickmode='array',
            tickvals=data['GENE'],
            ticktext=data['GENE'],
            automargin=True,
            tickfont=dict(size=10, color='black', family='Arial Black'),
            dtick=1
        )

        # Update x-axis with correct ratio range
        fig.update_xaxes(
            title_text='<b>Ratio</b>',
            title_font=dict(size=15, color='black', family='Arial Black'),
            row=1, col=i,
            range=[0, 1],
            tickvals=[0, 0.25, 0.5, 0.75, 1],
            ticktext=['0%', '25%', '50%', '75%', '100%'],
            tickfont=dict(size=12, color='black', family='Arial')
        )

    # Update overall layout
    fig.update_layout(
        title_text=f'<b>Heterozygosity Distribution for Chain {chain}</b>',
        title_font=dict(size=26, color='black', family='Arial Black'),
        title_x=0.5,
        height=subplot_height,
        width=1800,
        barmode='stack',
        bargap=0.15,
        margin=dict(l=150, r=50, t=130, b=50),
        uniformtext_minsize=8,
        uniformtext_mode='hide',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,
            xanchor="right",
            x=1,
            font=dict(size=13, color='black', family='Arial Black')
        ),
        plot_bgcolor='white',
    )

    # Adjust size of subplot titles
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=17, color='black', family='Arial Black')

    # Save as HTML
    fig.write_html(output_file, auto_open=False)