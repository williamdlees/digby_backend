import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from upsetplot import plot
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import json

__all__ = ['create_upset_plot']

def create_upset_plot(matrices, output_file, genes_to_plot=None):
    """
    Create UpSet plot for given matrices
    
    :param matrices: Dictionary of gene matrices
    :param output_file: Path to output file (.pdf or .html)
    :param genes_to_plot: Gene or list of genes to plot. If None, plots all genes. Defaults to None.
    :return: Boolean indicating success
    """
    file_ext = os.path.splitext(output_file)[1].lower()

    if file_ext not in ['.pdf', '.html']:
        raise ValueError("Output file must be either .pdf or .html")

    # Handle genes_to_plot filtering
    if genes_to_plot is not None:
        # Convert genes_to_plot to a list if it's a string
        if isinstance(genes_to_plot, str):
            genes_to_plot = [genes_to_plot]
        
        # Filter matrices to only include specified genes
        matrices = {k: v for k, v in matrices.items() if k in genes_to_plot}
        
        # Check if any matching genes were found
        if not matrices:
            print(f"No matching genes found. Available genes: {list(matrices.keys())}")
            return False

    try:
        if file_ext == '.pdf':
            return _create_pdf_upset(matrices, output_file)
        else:
            return _create_html_upset(matrices, output_file)
    except Exception as e:
        print(f"Error creating upset plot: {str(e)}")
        return False

    
def _create_pdf_upset(matrices, output_file):
    """Internal function to create PDF upset plots with sorted intersection sizes"""
    FIXED_WIDTH = 20
    FIXED_HEIGHT = 15
    SYMBOL_COLOR = '#0066CC'
    MAX_ALLELE_LENGTH = 40

    def clean_allele_name(name):
        """Clean allele names to ensure consistency"""
        # Remove any duplicate '01_' prefixes
        while '01_01_' in name:
            name = name.replace('01_01_', '01_')
        return name

    with PdfPages(output_file) as pdf:
        for gene, matrix in matrices.items():
            plt.clf()
            plt.close('all')

            fig = plt.figure()
            fig.set_size_inches(FIXED_WIDTH, FIXED_HEIGHT, forward=True)

            # Clean column names
            matrix.columns = [clean_allele_name(col) for col in matrix.columns]

            # Create symbols for long allele names
            allele_symbols = {}
            legend_needed = False
            for i, allele in enumerate(matrix.columns):
                if len(str(allele)) > MAX_ALLELE_LENGTH:
                    allele_symbols[allele] = f'A{i+1}'
                    legend_needed = True

            # Create a copy of matrix with symbolic names if needed
            if legend_needed:
                matrix_copy = matrix.copy()
                matrix_copy.columns = [allele_symbols.get(col, col) for col in matrix.columns]
                upset_data = matrix_copy.value_counts()
            else:
                upset_data = matrix.value_counts()

            # Sort the upset data by intersection sizes
            intersection_sizes = {}
            for pattern in upset_data.index:
                active_cols = [i for i, is_active in enumerate(pattern) if is_active]
                subset = matrix.iloc[:, active_cols]
                intersection_size = subset.all(axis=1).sum()
                intersection_sizes[pattern] = intersection_size

            # Sort the patterns by intersection size, and then by the pattern itself for consistency
            sorted_patterns = sorted(upset_data.index,
                                    key=lambda x: (-intersection_sizes[x], x))
            upset_data = upset_data.reindex(sorted_patterns)

            # Create the plot
            axes = plot(upset_data,
                        sort_by=None,
                        show_counts=True,
                        element_size=25 if len(matrix.columns) > 30 else 40,
                        fig=fig)

            if len(fig.axes) >= 2:
                intersection_ax = fig.axes[1]
                current_pos = intersection_ax.get_position()
                intersection_ax.set_position([0.8, current_pos.y0, current_pos.width, current_pos.height])

            fig.set_size_inches(FIXED_WIDTH, FIXED_HEIGHT, forward=True)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.85, bottom=0.15)

            if legend_needed:
                legend_text = []
                for allele, symbol in sorted(allele_symbols.items()):
                    legend_text.append(f'{symbol}: {allele}')

                plt.figtext(0.02, 0.02, 'Legend:\n' + '\n'.join(legend_text),
                            fontsize=7, family='monospace',
                            color=SYMBOL_COLOR,
                            bbox=dict(facecolor='white', alpha=0.8,
                                      edgecolor=SYMBOL_COLOR, linewidth=1),
                            verticalalignment='bottom')

            plt.suptitle(gene, y=0.95, fontsize=14, fontweight='bold')

            for ax in fig.axes:
                if ax.get_xaxis().get_visible():
                    ax.tick_params(axis='x', rotation=45, labelsize=10)
                ax.tick_params(axis='y', labelsize=10)

                for label in ax.get_yticklabels():
                    label_text = label.get_text()
                    if label_text.strip() in allele_symbols.values():
                        label.set_color(SYMBOL_COLOR)
                    label.set_fontsize(8)
                    label.set_fontweight('bold')

            pdf.savefig(fig,
                        dpi=300,
                        orientation='landscape',
                        bbox_inches=None)

            plt.close(fig)

    return True

# Previous imports and setup remain the same...

    # Previous imports and setup remain the same...
# Previous imports and setup remain the same...
# Previous imports and setup remain the same...
def _create_html_upset(matrices, output_file):
    """Create an HTML UpSet plot that matches the PDF output layout"""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pandas as pd
    import json
    
    html_start = '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>UpSet Plots</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            max-width: 2000px;
            margin: 0 auto; 
            padding: 20px;
            display: flex;
            flex-direction: column;
            align-items: center;
            background-color: #f8f8f8;
        }
        h1 {
            font-size: 36px;
            font-weight: bold;
            margin-bottom: 20px;
            color: #333;
        }
        #search-container {
            margin-bottom: 20px;
        }
        #plot-container {
            width: 100%;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 40px;
        }
        .gene-plot {
            width: 100%;
            max-width: 1800px;
            display: flex;
            justify-content: center;
            overflow-x: auto;
            background-color: white;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            padding: 20px;
        }
        #debug-info {
            background-color: #f0f0f0;
            padding: 10px;
            margin: 20px 0;
            white-space: pre-wrap;
            display: none;
            width: 100%;
            max-width: 1800px;
            border-radius: 5px;
            font-family: monospace;
            font-size: 14px;
        }
        #gene-search {
            padding: 10px;
            font-size: 16px;
            border: 1px solid #ccc;
            border-radius: 5px;
            margin-right: 10px;
        }
        button {
            padding: 10px 20px;
            font-size: 16px;
            background-color: #007bff;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
        }
        button:hover {
            background-color: #0056b3;
        }
    </style>
</head>
<body>
    <h1>UpSet Plots</h1>
    <div id="search-container">
        <label for="gene-search"><strong>Search Gene:</strong></label>
        <input type="text" id="gene-search" placeholder="Enter gene name">
        <button onclick="searchGene()">Search</button>
    </div>
    <div id="plot-container"></div>
    <div id="debug-info"></div>
    <script>
        // Debug logging function
        function debugLog(message) {
            var debugInfo = document.getElementById('debug-info');
            debugInfo.style.display = 'block';
            debugInfo.textContent += message + '\\n';
            console.log(message);
        }

        // Function to safely parse plot data
        function safeParsePlotData(plotlyData) {
            try {
                // Deep clone to avoid reference issues
                return JSON.parse(JSON.stringify({
                    data: plotlyData.data,
                    layout: plotlyData.layout
                }));
            } catch (err) {
                debugLog('Data parsing error: ' + err);
                return null;
            }
        }

        // Search function
        function searchGene() {
            var searchTerm = document.getElementById('gene-search').value.toLowerCase();
            var genePlots = document.getElementsByClassName('gene-plot');
            var firstMatchIndex = -1;

            for (var i = 0; i < genePlots.length; i++) {
                var geneName = genePlots[i].getAttribute('data-gene').toLowerCase();
                if (geneName.includes(searchTerm)) {
                    if (firstMatchIndex === -1) {
                        firstMatchIndex = i;
                    }
                    genePlots[i].style.display = 'flex';
                } else {
                    genePlots[i].style.display = 'flex';
                }
            }

            if (firstMatchIndex !== -1) {
                genePlots[firstMatchIndex].scrollIntoView({behavior: 'smooth', block: 'center'});
            }
        }

        // Main plot generation function
        function generateUpsetPlots(plotsData) {
            plotsData.forEach(function(plotData, index) {
                try {
                    var gene = plotData.layout.title.text.replace('<b>', '').replace('</b>', '');

                    // Create a container for each gene plot
                    var geneplotDiv = document.createElement('div');
                    geneplotDiv.className = 'gene-plot';
                    geneplotDiv.setAttribute('data-gene', gene);
                    
                    // Create a unique div for each plot
                    var plotDiv = document.createElement('div');
                    plotDiv.id = 'plot_' + index;
                    
                    geneplotDiv.appendChild(plotDiv);
                    document.getElementById('plot-container').appendChild(geneplotDiv);

                    // Safely parse plot data
                    var parsedData = safeParsePlotData(plotData);
                    if (!parsedData) {
                        throw new Error('Failed to parse plot data');
                    }

                    // Attempt to create plot
                    Plotly.newPlot(
                        plotDiv.id,
                        parsedData.data,
                        parsedData.layout,
                        {
                            responsive: true,
                            displayModeBar: true,
                            displaylogo: false,
                            modeBarButtonsToRemove: ['lasso2d', 'select2d']
                        }
                    ).catch(function(err) {
                        debugLog('Plotly plot error for ' + plotDiv.id + ': ' + err);
                        plotDiv.innerHTML = 'Error generating plot: ' + err;
                    });
                } catch (err) {
                    debugLog('Plot generation error: ' + err);
                }
            });
        }
    </script>
'''
    html_end = '''
</body>
</html>
'''
    html_content = [html_start]

    
    # Collect all plot data to pass to JavaScript
    plot_data_collection = []
    
    for gene, matrix in matrices.items():
        # Get unique alleles and sort them by frequency
        alleles = list(matrix.columns)
        allele_freqs = {allele: matrix[allele].sum() for allele in alleles}
        
        # Calculate pattern counts and intersection sizes
        patterns = []
        value_counts = matrix.value_counts()
        
        for pattern in value_counts.index:
            active_cols = [i for i, is_active in enumerate(pattern) if is_active]
            inactive_cols = [i for i, is_active in enumerate(pattern) if not is_active]
            
            if active_cols:
                pattern_count = value_counts[pattern]
                subset = matrix.iloc[:, active_cols]
                if inactive_cols:
                    inactive_subset = matrix.iloc[:, inactive_cols]
                    intersection_size = (subset.all(axis=1) & ~inactive_subset.any(axis=1)).sum()
                else:
                    intersection_size = subset.all(axis=1).sum()
                
                patterns.append({
                    'pattern': pattern,
                    'active_cols': active_cols,
                    'pattern_count': pattern_count,
                    'intersection_size': intersection_size
                })
        
        patterns.sort(key=lambda x: (-x['intersection_size'], -x['pattern_count']))
        
        # Adjust row height and plot layout more carefully
        max_alleles = len(alleles)
        row_height = max(1.5, 2 - (max_alleles * 0.05))  # Gradual reduction
        vertical_spacing = 0.1  # Consistent vertical spacing
        
        # Create figure with adjusted subplot configuration
        fig = make_subplots(
            rows=2, 
            cols=1,
            row_heights=[0.4, 0.6],  # Slight adjustment to proportions
            vertical_spacing=vertical_spacing,
            subplot_titles=('Set Sizes', 'Intersection Matrix')
        )
        
        # Add alternating background stripes
        for i in range(max_alleles):
            y_pos = (max_alleles - 1 - i) * row_height
            if i % 2 == 0:
                fig.add_shape(
                    type="rect",
                    xref="x2",
                    yref="y2",
                    x0=-0.5,
                    x1=len(patterns) - 0.5,
                    y0=y_pos - row_height/2,
                    y1=y_pos + row_height/2,
                    fillcolor="rgb(240,240,240)",
                    line_width=0,
                    layer="below"
                )
        
        # Add intersection size bars with text
        fig.add_trace(
            go.Bar(
                x=list(range(len(patterns))),
                y=[p['intersection_size'] for p in patterns],
                marker_color='rgb(53, 110, 196)',  # Richer blue
                showlegend=False,
                hovertemplate=(
                    "<b>Intersection Details</b><br>" +
                    "Size: %{y}<br>" +
                    "Alleles: %{customdata}<br>" +
                    "<extra></extra>"
                ),
                customdata=[", ".join([alleles[i] for i in p["active_cols"]]) for p in patterns],
                text=[str(p['intersection_size']) for p in patterns],
                textposition='outside',
                textfont=dict(
                    size=12,
                    color='rgb(30, 30, 30)',
                    family='Arial',
                    weight='bold'
                ),
                textangle=0,
                constraintext='none',
                cliponaxis=False
            ),
            row=1, col=1
        )
        
        # Draw dots for each allele row
        for i, allele in enumerate(alleles):
            y_pos = (max_alleles - 1 - i) * row_height
            
            dot_colors = []
            dot_sizes = []
            for pattern in patterns:
                col_idx = matrix.columns.get_loc(allele)
                is_active = col_idx in pattern['active_cols']
                dot_colors.append('black' if is_active else 'lightgray')
                dot_sizes.append(10 if is_active else 5)
            
            fig.add_trace(
                go.Scatter(
                    x=list(range(len(patterns))),
                    y=[y_pos] * len(patterns),
                    mode='markers',
                    marker=dict(
                        size=dot_sizes,
                        color=dot_colors,
                        line=dict(color='black', width=1)
                    ),
                    showlegend=False,
                    hoverinfo='text',
                    text=[f"{allele} (Total: {allele_freqs[allele]})" for _ in patterns]
                ),
                row=2, col=1
            )
        
        # Add vertical connections
        for x, pattern in enumerate(patterns):
            if len(pattern['active_cols']) > 1:
                y_positions = [(max_alleles - 1 - i) * row_height for i in pattern['active_cols']]
                fig.add_trace(
                    go.Scatter(
                        x=[x] * len(y_positions),
                        y=y_positions,
                        mode='lines',
                        line=dict(color='black', width=1),
                        showlegend=False,
                        hoverinfo='skip'
                    ),
                    row=2, col=1
                )
        
        # Update layout with more controlled sizing and increased margins
        fig.update_layout(
            title=dict(
                text=f"<b>{gene}</b>",
                x=0.5,
                y=0.99,
                xanchor='center',
                yanchor='top',
                font=dict(
                    size=24,
                    family='Arial Black',
                    color='rgb(30, 30, 30)'
                )
            ),
            showlegend=False,
            width=max(1200, min(1800, 100 * max(10, len(patterns)))),
            height=max(600, min(1000, 100 * max_alleles)),
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(l=400, r=100, t=150, b=100)  # Increased margins
        )
        
         # Update x-axes to ensure alignment
        fig.update_xaxes(matches='x', row=1, col=1)
        fig.update_xaxes(matches='x', row=2, col=1)
        fig.update_xaxes(showticklabels=False, showgrid=False, row=1, col=1)
        fig.update_xaxes(showticklabels=False, showgrid=False, row=2, col=1)
        
        # Update y-axes
        fig.update_yaxes(
            title=dict(
                text='<b>Intersection Size</b>',
                font=dict(size=14, family='Arial', color='rgb(30, 30, 30)')
            ),
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            row=1, col=1
        )
        
        y_labels = [f"{allele} ({allele_freqs[allele]})" for allele in alleles]
        y_tickvals = [i * row_height for i in range(max_alleles-1, -1, -1)]
        fig.update_yaxes(
            ticktext=y_labels,
            tickvals=y_tickvals,
            tickfont=dict(size=10, family='monospace'),
            showgrid=False,
            range=[-row_height, (max_alleles + 1) * row_height],
            row=2, col=1
        )
        
        # Add pattern count annotations
        annotations = []
        for i, allele in enumerate(alleles):
            y_pos = (max_alleles - 1 - i) * row_height
            count = sum(p['pattern_count'] for p in patterns 
                       if matrix.columns.get_loc(allele) in p['active_cols'])
            annotations.append(
                dict(
                    x=-0.15,
                    y=y_pos,
                    xref='paper',
                    yref='y2',
                    text=str(count),
                    showarrow=False,
                    font=dict(size=10),
                    align='right'
                )
            )
        
        fig.update_layout(annotations=annotations)
        
        # Collect plot data
        plot_data_collection.append({
            'data': fig.to_plotly_json()['data'],
            'layout': fig.to_plotly_json()['layout']
        })
    
    # Add script to generate plots
    html_content.append(f'''
    <script>
        // Generate plots when page loads
        document.addEventListener('DOMContentLoaded', function() {{
            try {{
                generateUpsetPlots({json.dumps(plot_data_collection)});
            }} catch (err) {{
                debugLog('Total plot generation error: ' + err);
            }}
        }});
    </script>
    ''')
    
    html_content.append(html_end)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(html_content))
    
    return True